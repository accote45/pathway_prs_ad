#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Compare PRSet competitive-pathway results for AD case/control between the
# APOE-EXCLUDED (original) run and the APOE-INCLUDED run.
#
# The two runs share every input except APOE handling (see
# scripts/apoe_included/README.md): the no-APOE run regresses APOE out of the
# phenotype, uses master_no_apoe.gmt, and excludes chr19:44409039-46412650;
# the with-APOE run keeps APOE in the phenotype, uses master.gmt, and clumps &
# scores the APOE region. The question here is the opposite of the subset
# analysis: we EXPECT the rankings to diverge, because the dominant chr19 APOE
# signal now flows into any pathway containing APOE-region genes. This produces:
#   1. logP concordance scatter - -log10(competitive P) no-APOE vs with-APOE,
#      permutation floor drawn on both axes, biggest movers + union top labelled.
#   2. R2 concordance scatter    - PRS R2 no-APOE vs with-APOE (effect size).
#   3. Top-gainers dumbbell       - pathways whose competitive signal increases
#      most when APOE is included (the APOE-driven pathways), no-APOE vs with.
#   4. Comparison CSV             - per-pathway P / R2 / rank in each run + deltas.
# Kendall tau-b / Spearman rho and top-N overlap are reported in the figure
# subtitles, the stats report, and the console.
#
# Deps: readr, dplyr, tidyr, stringr, ggplot2 (ggrepel used if available).
#
# Usage (from the repo root, with a tidyverse-capable R):
#   Rscript scripts/apoe_included/plot_prset_withapoe_vs_noapoe.R \
#     [--noapoe <no_apoe.summary>] [--withapoe <withapoe.summary>] \
#     [--outdir figures/withapoe_vs_noapoe] [--perm 10000] [--top 20]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2)
})

# ---- args (base-R, no optparse) -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) {
  i <- which(args == flag); if (length(i)) args[i + 1] else default
}
resultsdir <- getarg("--resultsdir",
  "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/results")
noapoe_path <- getarg("--noapoe",
  file.path(resultsdir, "ad_case.control_prset_nothreshold_eur.summary"))
withapoe_path <- getarg("--withapoe",
  file.path(resultsdir, "ad_case.control_prset_nothreshold_eur_withapoe.summary"))
outdir  <- getarg("--outdir", "figures/withapoe_vs_noapoe")
perm    <- as.integer(getarg("--perm", "10000"))
top     <- as.integer(getarg("--top", "20"))
godict  <- getarg("--godict",
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/gene_ontology/qced_data/aux_files/GO.dict")
mgidict <- getarg("--mgidict",
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/mgi/qced_data/MGI.dict")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
floor_p    <- 1 / (perm + 1)
floor_logp <- -log10(floor_p)
has_repel  <- requireNamespace("ggrepel", quietly = TRUE)

# ---- pretty pathway labels (GO/MGI dicts optional) ------------------------
load_dict <- function(path, label) {
  if (!file.exists(path)) {
    message(label, " dictionary not found (", path, "); those terms show IDs only.")
    return(character(0))
  }
  d <- read_tsv(path, show_col_types = FALSE, col_types = cols(.default = "c"))
  setNames(d$Name, d$ID)
}
term_map <- c(load_dict(godict, "GO"), load_dict(mgidict, "MGI"))

prettify_set <- function(set) {
  parts <- str_split_fixed(set, "__", 2)
  db <- parts[, 1]; name <- parts[, 2]
  no_sep <- name == ""
  db[no_sep] <- ""; name[no_sep] <- set[no_sep]
  pure_id <- str_detect(name, "^(GO|MP):[0-9]+$")
  hit <- pure_id & name %in% names(term_map)
  name[hit] <- unname(term_map[name[hit]])
  name <- str_replace_all(name, "_", " ")
  msig <- db == "MSigDB"
  name[msig] <- str_to_title(name[msig])
  ifelse(db == "", name, paste0(db, ": ", name))
}

# ---- read + harmonise -----------------------------------------------------
read_prset <- function(path, run_label) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read_table(path, show_col_types = FALSE)
  r2_col <- if ("PRS.R2" %in% names(df)) "PRS.R2"
            else if ("R2" %in% names(df)) "R2" else NA_character_
  if (is.na(r2_col)) stop("No R2 column in ", path)
  if (!"Competitive.P" %in% names(df))
    stop("No Competitive.P column in ", path, " (was --set-perm used?)")
  df %>%
    rename(R2 = all_of(r2_col)) %>%
    filter(!Set %in% c("Base", "Background")) %>%
    mutate(Run = run_label,
           logP       = -log10(pmax(Competitive.P, floor_p)),
           at_floor   = Competitive.P <= floor_p + 1e-12,
           R2_per_SNP = R2 / Num_SNP) %>%
    # consistent ranking rule across runs: comp P asc, then per-SNP density
    arrange(Competitive.P, desc(R2_per_SNP)) %>%
    mutate(rank = row_number())
}

noapoe   <- read_prset(noapoe_path,   "APOE excluded")
withapoe <- read_prset(withapoe_path, "APOE included")

# per-pathway wide table (inner join = pathways scored in both runs)
key_cols <- c("Set", "Competitive.P", "logP", "at_floor", "R2", "Coefficient",
              "Num_SNP", "rank")
wide <- inner_join(
  noapoe   %>% select(all_of(key_cols)) %>% rename_with(~paste0(., "_noapoe"),   -Set),
  withapoe %>% select(all_of(key_cols)) %>% rename_with(~paste0(., "_withapoe"), -Set),
  by = "Set"
) %>% mutate(Set_pretty  = prettify_set(Set),
             logP_delta  = logP_withapoe - logP_noapoe,   # + = gains signal w/ APOE
             rank_delta  = rank_withapoe - rank_noapoe,    # - = moves up w/ APOE
             SNP_delta   = Num_SNP_withapoe - Num_SNP_noapoe)

n_noapoe_only   <- nrow(noapoe)   - nrow(wide)
n_withapoe_only <- nrow(withapoe) - nrow(wide)
if (n_noapoe_only + n_withapoe_only > 0)
  message(sprintf("Note: %d pathway(s) only in no-APOE, %d only in with-APOE; ",
                  n_noapoe_only, n_withapoe_only),
          "concordance figures use the ", nrow(wide), " scored in both.")

# ---- summary statistics ---------------------------------------------------
top_noapoe   <- noapoe   %>% slice_head(n = top) %>% pull(Set)
top_withapoe <- withapoe %>% slice_head(n = top) %>% pull(Set)
overlap      <- length(intersect(top_noapoe, top_withapoe))

rho_logP <- suppressWarnings(cor(wide$logP_noapoe, wide$logP_withapoe,
                                 method = "spearman", use = "complete.obs"))
rho_R2   <- suppressWarnings(cor(wide$R2_noapoe, wide$R2_withapoe,
                                 method = "spearman", use = "complete.obs"))
rho_rank <- suppressWarnings(cor(wide$rank_noapoe, wide$rank_withapoe,
                                 method = "spearman", use = "complete.obs"))

# Kendall tau-b (tie-aware) over ALL pathways scored in both runs.
tau_logP <- suppressWarnings(cor(wide$logP_noapoe, wide$logP_withapoe,
                                 method = "kendall", use = "complete.obs"))
tau_R2   <- suppressWarnings(cor(wide$R2_noapoe, wide$R2_withapoe,
                                 method = "kendall", use = "complete.obs"))

# Top-K replication: overlap of each run's top K + one-sided hypergeometric P.
N_both <- nrow(wide)
overlap_test <- function(K) {
  K  <- min(K, N_both)
  tn <- wide %>% arrange(rank_noapoe)   %>% slice_head(n = K) %>% pull(Set)
  tw <- wide %>% arrange(rank_withapoe) %>% slice_head(n = K) %>% pull(Set)
  ov <- length(intersect(tn, tw))
  p  <- phyper(ov - 1, m = K, n = N_both - K, k = K, lower.tail = FALSE)
  data.frame(K = K, overlap = ov, expected = K * K / N_both,
             jaccard = ov / (2 * K - ov), hyper_p = p)
}
topK_tbl <- do.call(rbind, lapply(c(50, 100, 200), overlap_test))
h100 <- topK_tbl[which.min(abs(topK_tbl$K - 100)), ]

wide <- wide %>% mutate(is_top = Set %in% union(top_noapoe, top_withapoe))

# Biggest movers: pathways gaining the most competitive signal with APOE. These
# are the pathways where the chr19 APOE region enters and lifts the score; label
# them (plus the union of each run's top-N) on the concordance scatter.
top_gainers <- wide %>% arrange(desc(logP_delta)) %>% slice_head(n = top) %>% pull(Set)
wide <- wide %>% mutate(is_gainer = Set %in% top_gainers)

# Scatter labels: the gainers pile onto the permutation floor (y = ceiling), so
# labelling all of them is unreadable. Label only a small set - each run's top
# few pathways - and let the dumbbell (Fig 3) + CSV carry the full gainer list.
n_lab       <- 6L
lab_sets    <- union(head(top_noapoe, n_lab), head(top_withapoe, n_lab))
wide        <- wide %>% mutate(label_me = Set %in% lab_sets)

# reusable labeller
add_labels <- function(p, data, x, y, which_col = "label_me") {
  lab <- data[data[[which_col]], ]
  if (has_repel) {
    p + ggrepel::geom_text_repel(data = lab, aes({{x}}, {{y}}, label = Set_pretty),
          size = 2.4, max.overlaps = 25, colour = "grey20", min.segment.length = 0,
          box.padding = 0.5, seed = 1)
  } else {
    p + geom_text(data = lab, aes({{x}}, {{y}}, label = Set_pretty),
          size = 2.4, colour = "grey20", vjust = -0.6, check_overlap = TRUE)
  }
}

# ---- Figure 1: -log10(competitive P) concordance --------------------------
lim <- c(0, max(wide$logP_noapoe, wide$logP_withapoe) * 1.05)
p1 <- ggplot(wide, aes(logP_noapoe, logP_withapoe)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = floor_logp, linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = floor_logp, linetype = "dotted", colour = "grey60") +
  geom_point(aes(colour = is_gainer, alpha = is_gainer)) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "top APOE gainer", `FALSE` = "other"),
                      name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression(-log[10]~"(competitive P), APOE excluded"),
       y = expression(-log[10]~"(competitive P), APOE included"),
       title = "PRSet pathway signal: APOE included vs excluded (AD case/control)",
       subtitle = sprintf("Kendall tau-b (all %d paths) = %.2f; top-100 overlap = %d (exp %.1f), hypergeometric P = %.1e",
                          N_both, tau_logP, h100$overlap, h100$expected, h100$hyper_p),
       caption = "Points above the diagonal gain signal when APOE is included. Dotted lines = permutation floor (1/(perm+1)).") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")
p1 <- add_labels(p1, wide, logP_noapoe, logP_withapoe)
ggsave(file.path(outdir, "prset_withapoe_vs_noapoe_logP.pdf"), p1, width = 7.5, height = 7.5)

# ---- Figure 2: PRS R2 concordance -----------------------------------------
p2 <- ggplot(wide, aes(R2_noapoe, R2_withapoe)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = is_gainer, alpha = is_gainer)) +
  scale_colour_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "top APOE gainer", `FALSE` = "other"),
                      name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  labs(x = expression(PRS~R^2~", APOE excluded"), y = expression(PRS~R^2~", APOE included"),
       title = "PRSet pathway effect sizes: APOE included vs excluded",
       subtitle = sprintf("Spearman rho (R2) = %.2f over %d pathways scored in both runs",
                          rho_R2, N_both)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")
p2 <- add_labels(p2, wide, R2_noapoe, R2_withapoe)
ggsave(file.path(outdir, "prset_withapoe_vs_noapoe_R2.pdf"), p2, width = 7.5, height = 7.5)

# ---- Figure 3: top-gainers dumbbell ---------------------------------------
# Pathways whose competitive signal increases most when APOE is included.
dumb <- wide %>%
  filter(Set %in% top_gainers) %>%
  arrange(logP_delta) %>%
  mutate(Set_pretty = factor(Set_pretty, levels = Set_pretty)) %>%
  select(Set_pretty,
         `APOE excluded` = logP_noapoe, `APOE included` = logP_withapoe) %>%
  pivot_longer(c(`APOE excluded`, `APOE included`),
               names_to = "Run", values_to = "logP")

p3 <- ggplot(dumb, aes(logP, Set_pretty)) +
  geom_vline(xintercept = floor_logp, linetype = "dotted", colour = "grey60") +
  geom_line(aes(group = Set_pretty), colour = "grey75", linewidth = 0.6) +
  geom_point(aes(colour = Run), size = 2.6) +
  scale_colour_manual(values = c(`APOE excluded` = "grey55", `APOE included` = "#D55E00")) +
  labs(x = expression(-log[10]~"(competitive P)"), y = NULL,
       title = sprintf("Top %d pathways gaining signal when APOE is included", top),
       subtitle = "Dotted line = permutation floor",
       colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), legend.position = "top")
ggsave(file.path(outdir, "prset_withapoe_vs_noapoe_gainers_dumbbell.pdf"), p3,
       width = 9, height = max(4, 0.32 * top + 2), limitsize = FALSE)

# ---- comparison table -----------------------------------------------------
wide %>%
  transmute(Set, Set_pretty,
            CompP_noapoe = Competitive.P_noapoe, CompP_withapoe = Competitive.P_withapoe,
            logP_noapoe, logP_withapoe, logP_delta,
            R2_noapoe, R2_withapoe,
            Coef_noapoe = Coefficient_noapoe, Coef_withapoe = Coefficient_withapoe,
            Num_SNP_noapoe, Num_SNP_withapoe, SNP_delta,
            rank_noapoe, rank_withapoe, rank_delta,
            top_noapoe = Set %in% top_noapoe, top_withapoe = Set %in% top_withapoe,
            top_gainer = Set %in% top_gainers) %>%
  arrange(desc(logP_delta)) %>%
  write_csv(file.path(outdir, "prset_withapoe_vs_noapoe_table.csv"))

# ---- statistics report ----------------------------------------------------
stats_path <- file.path(outdir, "prset_withapoe_vs_noapoe_stats.txt")
writeLines(c(
  "PRSet APOE-included vs APOE-excluded concordance statistics (AD case/control)",
  "============================================================================",
  sprintf("APOE excluded : %s (%d pathways)", noapoe_path, nrow(noapoe)),
  sprintf("APOE included : %s (%d pathways)", withapoe_path, nrow(withapoe)),
  sprintf("Pathways scored in both (used below): %d", N_both),
  "",
  "1. Global rank concordance over ALL pathways",
  sprintf("   Kendall tau-b : logP = %+.3f | R2 = %+.3f", tau_logP, tau_R2),
  sprintf("   Spearman rho  : logP = %+.3f | R2 = %+.3f | rank = %+.3f",
          rho_logP, rho_R2, rho_rank),
  "",
  "2. Top-K replication: overlap of each run's top K + one-sided hypergeometric P",
  sprintf("   K=%-4d overlap=%-4d expected=%6.2f jaccard=%.3f hyper_P=%.3e",
          topK_tbl$K, topK_tbl$overlap, topK_tbl$expected,
          topK_tbl$jaccard, topK_tbl$hyper_p),
  "",
  sprintf("3. Top-%d union overlap between runs: %d pathway(s)", top, overlap),
  "",
  "Notes",
  " - Unlike the subset analysis, LOW concordance is expected and informative:",
  "   including APOE routes the dominant chr19 signal into every pathway that",
  "   contains an APOE-region gene, so those pathways gain competitive signal.",
  " - Figure 3 / the CSV (sorted by logP_delta) list the pathways that gain the",
  "   most; these are the APOE-driven sets to interpret with caution.",
  " - Competitive P is floored at 1/(perm+1); heavy ties at the floor limit any",
  "   P-based correlation. Kendall tau-b is tie-aware; the top-K hypergeometric",
  "   test is floor-robust (set membership, not fine ranks)."
), stats_path)
message("Wrote stats report to ", stats_path)

message("Done. Wrote 3 figures + prset_withapoe_vs_noapoe_table.csv + stats report to ", outdir)
message(sprintf("Kendall tau-b (logP) = %.3f | Spearman rho (logP) = %.3f | top-100 overlap = %d (exp %.1f, P = %.2e)",
                tau_logP, rho_logP, h100$overlap, h100$expected, h100$hyper_p))
if (!has_repel)
  message("ggrepel not installed: point labels use geom_text (may overlap). ",
          "install.packages('ggrepel') for cleaner labels.")
