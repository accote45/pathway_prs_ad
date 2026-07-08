#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Plot PRSet competitive-pathway results for AD age-of-onset and AD diagnosis.
#
# Handles the permutation-floor problem: with --set-perm N the smallest
# competitive P is 1/(N+1), so many top pathways tie at the ceiling of
# -log10(Competitive.P). Rather than hide that, the main figure draws the
# floor explicitly and uses PRS.R2 (colour) + Num_SNP (size) to reveal the
# real ordering among tied pathways.
#
# Usage:
#   Rscript plot_prset_results.R \
#     --aoo  path/to/ad_ageofonset_prset_..._eur.summary \
#     --addx path/to/ad_dx_prset_..._eur.summary \
#     --perm 10000 --top 20 --outdir figures/
#
# All flags are optional; defaults point at the cluster results paths.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ---- args -----------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = list(
  make_option("--aoo",    type = "character",
              default = "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/results/ad_ageofonset_prset_nothreshold_eur.summary",
              help = "PRSet .summary file for AD age of onset"),
  make_option("--addx",   type = "character",
              default = "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/results/ad_case.control_prset_nothreshold_eur.summary",
              help = "PRSet .summary file for AD diagnosis"),
  make_option("--perm",   type = "integer",   default = 10000L,
              help = "--set-perm value used in PRSet [default %default]"),
  make_option("--top",    type = "integer",   default = 15L,
              help = "Number of top pathways per phenotype [default %default]"),
  make_option("--godict", type = "character",
              default = "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/gene_ontology/qced_data/aux_files/GO.dict",
              help = "GO ID -> term name dictionary (tab-separated: ID, Name)"),
  make_option("--mgidict", type = "character",
              default = "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/mgi/qced_data/MGI.dict",
              help = "MGI (MP) ID -> term name dictionary (tab-separated: ID, Name)"),
  make_option("--outdir", type = "character", default = "figures",
              help = "Output directory [default %default]")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
floor_p    <- 1 / (opt$perm + 1)          # smallest achievable competitive P
floor_logp <- -log10(floor_p)

# ---- pretty pathway labels ------------------------------------------------
# ID -> term name lookup (optional; falls back to the raw ID if missing).
# Both dicts are tab-separated with columns ID, Name; GO uses "GO:####" IDs
# and MGI uses "MP:####" IDs, so a single combined map is unambiguous.
load_dict <- function(path, label) {
  if (!file.exists(path)) {
    message(label, " dictionary not found (", path,
            "); those terms will show IDs only.")
    return(character(0))
  }
  d <- readr::read_tsv(path, show_col_types = FALSE,
                       col_types = readr::cols(.default = "c"))
  setNames(d$Name, d$ID)
}
term_map <- c(load_dict(opt$godict, "GO"), load_dict(opt$mgidict, "MGI"))

# Turn "MSigDB__BORCZUK_MALIGNANT_MESOTHELIOMA_UP" into
# "MSigDB: BORCZUK MALIGNANT MESOTHELIOMA UP", a bare GO term
# "GO__GO:0043032" into "GO: mitochondrion inheritance", and an MGI term
# "MGI__MP:0000547" into "MGI: <phenotype description>".
prettify_set <- function(set) {
  parts <- str_split_fixed(set, "__", 2)
  db    <- parts[, 1]
  name  <- parts[, 2]
  no_sep <- name == ""                 # no "__" separator present
  db[no_sep]   <- ""
  name[no_sep] <- set[no_sep]
  # Pure ID terms ("GO:0043032" / "MP:0000547"): swap the ID for its name.
  pure_id <- str_detect(name, "^(GO|MP):[0-9]+$")
  hit     <- pure_id & name %in% names(term_map)
  name[hit] <- unname(term_map[name[hit]])
  name <- str_replace_all(name, "_", " ")
  # MSigDB gene-set names are ALL CAPS; title-case them for readability.
  msig <- db == "MSigDB"
  name[msig] <- str_to_title(name[msig])
  ifelse(db == "", name, paste0(db, ": ", name))
}

# ---- read + harmonise -----------------------------------------------------
# PRSet .summary columns vary slightly by version; detect the R2 column.
read_prset <- function(path, pheno_label) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- readr::read_table(path, show_col_types = FALSE)

  r2_col <- dplyr::case_when(
    "PRS.R2" %in% names(df) ~ "PRS.R2",
    "R2"     %in% names(df) ~ "R2",
    TRUE ~ NA_character_
  )
  if (is.na(r2_col)) stop("No R2 column found in ", path)
  if (!"Competitive.P" %in% names(df))
    stop("No Competitive.P column in ", path, " (was --set-perm used?)")

  df %>%
    rename(R2 = all_of(r2_col)) %>%
    mutate(Phenotype = pheno_label) %>%
    # drop the genome-wide background set; keep only pathways
    filter(!Set %in% c("Base", "Background"))
}

dat <- bind_rows(
  read_prset(opt$aoo,  "AD age of onset"),
  read_prset(opt$addx, "AD diagnosis")
)

# ---- composite ranking ----------------------------------------------------
# Competitive.P (asc) -> R2 / Num_SNP (desc, per-SNP density) as tiebreaker.
ranked <- dat %>%
  mutate(
    at_floor    = Competitive.P <= floor_p + 1e-12,
    logP        = -log10(pmax(Competitive.P, floor_p)),  # cap at the ceiling
    R2_per_SNP  = R2 / Num_SNP
  ) %>%
  group_by(Phenotype) %>%
  arrange(Competitive.P, desc(R2_per_SNP), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  slice_head(n = opt$top) %>%
  ungroup()

# per-facet ordering (reorder_within, inlined to avoid a tidytext dependency)
ranked <- ranked %>%
  mutate(Set_ord = paste(Set, Phenotype, sep = "___"))
lvls <- ranked %>% arrange(Phenotype, desc(rank)) %>% pull(Set_ord)
ranked$Set_ord <- factor(ranked$Set_ord, levels = unique(lvls))
strip_ord <- function(x) prettify_set(str_remove(x, "___.*$"))

# ---- Figure 1: faceted lollipop ------------------------------------------
p1 <- ggplot(ranked, aes(x = logP, y = Set_ord)) +
  geom_vline(xintercept = floor_logp, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = logP, yend = Set_ord),
               colour = "grey80", linewidth = 0.5) +
  geom_point(aes(colour = R2, size = Num_SNP)) +
  scale_y_discrete(labels = strip_ord) +
  scale_colour_viridis_c(option = "C", name = expression(PRS~R^2)) +
  scale_size_continuous(name = "N SNPs", range = c(2, 7)) +
  facet_wrap(~ Phenotype, scales = "free_y") +
  labs(x = expression(-log[10]~"(competitive P)"), y = NULL,
       title = sprintf("Top %d pathways by competitive P", opt$top)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(opt$outdir, "prset_lollipop_faceted.pdf"),
       p1, width = 11, height = 6)

# ---- Figure 2: cross-phenotype concordance --------------------------------
# R2 for every pathway in one phenotype vs the other; label the union of tops.
top_sets <- unique(ranked$Set)
conc <- dat %>%
  select(Phenotype, Set, R2, Competitive.P) %>%
  pivot_wider(names_from = Phenotype,
              values_from = c(R2, Competitive.P))

aoo_r2  <- grep("^R2_AD age", names(conc), value = TRUE)
addx_r2 <- grep("^R2_AD diag", names(conc), value = TRUE)

if (length(aoo_r2) == 1 && length(addx_r2) == 1) {
  conc <- conc %>%
    rename(R2_aoo = all_of(aoo_r2), R2_addx = all_of(addx_r2)) %>%
    mutate(is_top = Set %in% top_sets,
           Set_label = prettify_set(Set))

  p2 <- ggplot(conc, aes(R2_aoo, R2_addx)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(aes(colour = is_top, alpha = is_top)) +
    ggrepel::geom_text_repel(
      data = subset(conc, is_top), aes(label = Set_label),
      size = 2.6, max.overlaps = 20, colour = "grey20") +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70"),
                        guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
    labs(x = expression(PRS~R^2~" (AD age of onset)"),
         y = expression(PRS~R^2~" (AD diagnosis)"),
         title = "Cross-phenotype concordance of pathway PRS effects") +
    theme_minimal(base_size = 11)

  ggsave(file.path(opt$outdir, "prset_concordance.pdf"),
         p2, width = 6.5, height = 6)
} else {
  message("Skipped concordance plot: could not resolve both R2 columns.")
}

# ---- ranked table (also useful as a supplement) ---------------------------
ranked %>%
  transmute(Phenotype, rank, Set, Competitive.P, at_floor,
            R2, R2_per_SNP, Coefficient, Num_SNP) %>%
  arrange(Phenotype, rank) %>%
  write_csv(file.path(opt$outdir, "prset_top_pathways.csv"))

message("Done. Wrote figures + prset_top_pathways.csv to ", opt$outdir)
if (any(ranked$at_floor))
  message(sprintf("Note: %d of the top pathways are tied at the permutation floor. ",
                  sum(ranked$at_floor)),
          "Consider re-running the top tier with a larger --set-perm to break ties.")
