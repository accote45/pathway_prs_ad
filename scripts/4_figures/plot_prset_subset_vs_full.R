#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Compare PRSet competitive-pathway results for AD case/control between the
# FULL EUR sample and the 1000-case/1000-control SUBSET run.
#
# Both runs use identical inputs (base GWAS, master_no_apoe.gmt, background,
# chr19 APOE exclusion, --set-perm 10000); the subset differs only in --keep.
# The question is whether the down-sampled run recovers the full-sample
# pathway signal. This produces:
#   1. logP concordance scatter  - -log10(competitive P) full vs subset, with
#      the permutation floor drawn on both axes and top pathways labelled.
#   2. R2 concordance scatter    - PRS R2 full vs subset (no floor; effect-size
#      replication).
#   3. Top-pathway dumbbell      - for the full run's top-N pathways, full vs
#      subset -log10(competitive P) connected, to see which lose signal.
#   4. Comparison CSV            - per-pathway P / R2 / rank in each run + deltas.
# Spearman correlations (P-rank, R2) and top-N overlap are reported in the
# figure subtitles and the console.
#
# Deps: readr, dplyr, tidyr, stringr, ggplot2 (ggrepel used if available).
#
# Usage (on Minerva, from the repo root, with a tidyverse-capable R):
#   Rscript scripts/plot_prset_subset_vs_full.R \
#     [--full <full.summary>] [--subset <subset.summary>] \
#     [--outdir figures/subset_vs_full] [--perm 10000] [--top 15]
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
full_path <- getarg("--full",
  file.path(resultsdir, "ad_case.control_prset_nothreshold_eur.summary"))
subset_path <- getarg("--subset",
  file.path(resultsdir, "ad_case.control_prset_nothreshold_eur_subset1000.summary"))
outdir  <- getarg("--outdir", "figures/subset_vs_full")
perm    <- as.integer(getarg("--perm", "10000"))
top     <- as.integer(getarg("--top", "15"))
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

full   <- read_prset(full_path,   "Full")
subset <- read_prset(subset_path, "Subset (1000/1000)")

# per-pathway wide table (inner join = pathways scored in both runs)
key_cols <- c("Set", "Competitive.P", "logP", "at_floor", "R2", "Coefficient",
              "Num_SNP", "rank")
wide <- inner_join(
  full   %>% select(all_of(key_cols)) %>% rename_with(~paste0(., "_full"),   -Set),
  subset %>% select(all_of(key_cols)) %>% rename_with(~paste0(., "_subset"), -Set),
  by = "Set"
) %>% mutate(Set_pretty = prettify_set(Set))

n_full_only   <- nrow(full)   - nrow(wide)
n_subset_only <- nrow(subset) - nrow(wide)
if (n_full_only + n_subset_only > 0)
  message(sprintf("Note: %d pathway(s) only in full, %d only in subset; ",
                  n_full_only, n_subset_only),
          "concordance figures use the ", nrow(wide), " scored in both.")

# ---- summary statistics ---------------------------------------------------
top_full   <- full   %>% slice_head(n = top) %>% pull(Set)
top_subset <- subset %>% slice_head(n = top) %>% pull(Set)
overlap    <- length(intersect(top_full, top_subset))

rho_logP <- suppressWarnings(cor(wide$logP_full, wide$logP_subset,
                                 method = "spearman", use = "complete.obs"))
rho_R2   <- suppressWarnings(cor(wide$R2_full, wide$R2_subset,
                                 method = "spearman", use = "complete.obs"))
rho_rank <- suppressWarnings(cor(wide$rank_full, wide$rank_subset,
                                 method = "spearman", use = "complete.obs"))

# Kendall tau-b (tie-aware) over ALL pathways scored in both runs. This is the
# global concordance statistic; run over all pathways (not a selected top-K) to
# avoid range-restriction bias.
tau_logP <- suppressWarnings(cor(wide$logP_full, wide$logP_subset,
                                 method = "kendall", use = "complete.obs"))
tau_R2   <- suppressWarnings(cor(wide$R2_full, wide$R2_subset,
                                 method = "kendall", use = "complete.obs"))

# Top-K replication: overlap of each run's top K + one-sided hypergeometric
# (Fisher) test. Floor-robust because it uses set membership, not fine ranks.
N_both <- nrow(wide)
overlap_test <- function(K) {
  K  <- min(K, N_both)
  tf <- wide %>% arrange(rank_full)   %>% slice_head(n = K) %>% pull(Set)
  ts <- wide %>% arrange(rank_subset) %>% slice_head(n = K) %>% pull(Set)
  ov <- length(intersect(tf, ts))
  # P(overlap >= ov): white = full top-K (K), black = N - K, draws = subset top-K (K)
  p  <- phyper(ov - 1, m = K, n = N_both - K, k = K, lower.tail = FALSE)
  data.frame(K = K, overlap = ov, expected = K * K / N_both,
             jaccard = ov / (2 * K - ov), hyper_p = p)
}
topK_tbl <- do.call(rbind, lapply(c(50, 100, 200), overlap_test))
h100 <- topK_tbl[which.min(abs(topK_tbl$K - 100)), ]

# Rank-enrichment (the correct KS-style test): are the full run's top-K pathways
# positioned higher in the SUBSET ranking than the rest? AUC = P(a full-top
# pathway outranks a random non-top pathway in the subset); 0.5 = chance.
K_enrich   <- min(100, N_both)
top_full_K <- wide %>% arrange(rank_full) %>% slice_head(n = K_enrich) %>% pull(Set)
enr_grp    <- wide$Set %in% top_full_K
wt <- suppressWarnings(wilcox.test(wide$logP_subset[enr_grp],
                                   wide$logP_subset[!enr_grp],
                                   alternative = "greater"))
auc_enrich <- unname(wt$statistic) / (sum(enr_grp) * sum(!enr_grp))
p_enrich   <- wt$p.value

wide <- wide %>% mutate(is_top = Set %in% union(top_full, top_subset))

# reusable labeller for the top pathways
add_labels <- function(p, data, x, y) {
  lab <- subset(data, is_top)
  if (has_repel) {
    p + ggrepel::geom_text_repel(data = lab, aes({{x}}, {{y}}, label = Set_pretty),
          size = 2.4, max.overlaps = 20, colour = "grey20", min.segment.length = 0)
  } else {
    p + geom_text(data = lab, aes({{x}}, {{y}}, label = Set_pretty),
          size = 2.4, colour = "grey20", vjust = -0.6, check_overlap = TRUE)
  }
}

# ---- Figure 1: -log10(competitive P) concordance --------------------------
lim <- c(0, max(wide$logP_full, wide$logP_subset) * 1.05)
p1 <- ggplot(wide, aes(logP_full, logP_subset)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = floor_logp, linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = floor_logp, linetype = "dotted", colour = "grey60") +
  geom_point(aes(colour = is_top, alpha = is_top)) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70"), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression(-log[10]~"(competitive P), full"),
       y = expression(-log[10]~"(competitive P), subset"),
       title = "PRSet pathway signal: subset vs full (AD case/control)",
       subtitle = sprintf("Kendall tau-b (all %d paths) = %.2f; top-100 overlap = %d (exp %.1f), hypergeometric P = %.1e",
                          N_both, tau_logP, h100$overlap, h100$expected, h100$hyper_p),
       caption = "Dotted lines = permutation floor (1/(perm+1)); labelled = union of each run's top pathways") +
  theme_minimal(base_size = 11)
p1 <- add_labels(p1, wide, logP_full, logP_subset)
ggsave(file.path(outdir, "prset_subset_vs_full_logP.pdf"), p1, width = 7.5, height = 7)

# ---- Figure 2: PRS R2 concordance (no floor) ------------------------------
p2 <- ggplot(wide, aes(R2_full, R2_subset)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = is_top, alpha = is_top)) +
  scale_colour_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "grey70"), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
  labs(x = expression(PRS~R^2~", full"), y = expression(PRS~R^2~", subset"),
       title = "PRSet pathway effect sizes: subset vs full",
       subtitle = sprintf("Spearman rho = %.2f (subset R2 is inflated by small N; rank concordance is the point)",
                          rho_R2)) +
  theme_minimal(base_size = 11)
p2 <- add_labels(p2, wide, R2_full, R2_subset)
ggsave(file.path(outdir, "prset_subset_vs_full_R2.pdf"), p2, width = 7.5, height = 7)

# ---- Figure 3: top-pathway dumbbell ---------------------------------------
dumb <- wide %>%
  filter(Set %in% top_full) %>%
  arrange(logP_full) %>%
  mutate(Set_pretty = factor(Set_pretty, levels = Set_pretty)) %>%
  select(Set_pretty, Full = logP_full, `Subset (1000/1000)` = logP_subset) %>%
  pivot_longer(c(Full, `Subset (1000/1000)`), names_to = "Run", values_to = "logP")

p3 <- ggplot(dumb, aes(logP, Set_pretty)) +
  geom_vline(xintercept = floor_logp, linetype = "dotted", colour = "grey60") +
  geom_line(aes(group = Set_pretty), colour = "grey75", linewidth = 0.6) +
  geom_point(aes(colour = Run), size = 2.6) +
  scale_colour_manual(values = c(`Full` = "#D55E00", `Subset (1000/1000)` = "#0072B2")) +
  labs(x = expression(-log[10]~"(competitive P)"), y = NULL,
       title = sprintf("Full run's top %d pathways: signal in the subset", top),
       subtitle = "Dotted line = permutation floor",
       colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), legend.position = "top")
ggsave(file.path(outdir, "prset_subset_vs_full_dumbbell.pdf"), p3,
       width = 9, height = max(4, 0.32 * top + 2), limitsize = FALSE)

# ---- comparison table -----------------------------------------------------
wide %>%
  transmute(Set, Set_pretty,
            CompP_full = Competitive.P_full, CompP_subset = Competitive.P_subset,
            logP_full, logP_subset,
            R2_full, R2_subset, Coef_full = Coefficient_full, Coef_subset = Coefficient_subset,
            Num_SNP_full, Num_SNP_subset,
            rank_full, rank_subset, rank_delta = rank_subset - rank_full,
            top_full = Set %in% top_full, top_subset = Set %in% top_subset) %>%
  arrange(rank_full) %>%
  write_csv(file.path(outdir, "prset_subset_vs_full_table.csv"))

# ---- statistics report ----------------------------------------------------
stats_path <- file.path(outdir, "prset_subset_vs_full_stats.txt")
writeLines(c(
  "PRSet subset-1000 vs full concordance statistics (AD case/control)",
  "===================================================================",
  sprintf("Full   : %s (%d pathways)", full_path, nrow(full)),
  sprintf("Subset : %s (%d pathways)", subset_path, nrow(subset)),
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
  sprintf("3. Rank enrichment (Mann-Whitney): full run's top-%d in the subset ranking",
          K_enrich),
  sprintf("   AUC = %.3f (0.5 = chance), one-sided P = %.3e", auc_enrich, p_enrich),
  "",
  "Notes",
  " - Competitive P is floored at 1/(perm+1); heavy ties at the floor limit any",
  "   P-based correlation. Kendall tau-b is tie-aware; the top-K hypergeometric",
  "   test is floor-robust (set membership, not fine ranks) -> lead with it.",
  " - R2 is NOT gene-set-size adjusted (competitive P is); R2 and P concordance",
  "   answer subtly different questions.",
  " - A two-sample KS on the marginal P/R2 distributions was deliberately NOT",
  "   used: it is invariant to pathway labels and tests distribution shape, not",
  "   concordance."
), stats_path)
message("Wrote stats report to ", stats_path)

message("Done. Wrote 3 figures + prset_subset_vs_full_table.csv + stats report to ", outdir)
message(sprintf("Kendall tau-b (logP) = %.3f | Spearman rho (logP) = %.3f | top-100 overlap = %d (exp %.1f, P = %.2e) | enrichment AUC = %.3f",
                tau_logP, rho_logP, h100$overlap, h100$expected, h100$hyper_p, auc_enrich))
if (!has_repel)
  message("ggrepel not installed: point labels use geom_text (may overlap). ",
          "install.packages('ggrepel') for cleaner labels.")
