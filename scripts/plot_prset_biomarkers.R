#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Plot PRSet competitive-pathway results across the AD biomarker panel:
# plasma proteins (GFAP, NEFL) and MRI measures (hippocampal volume, WMH).
#
# Generalizes plot_prset_results.R (which is hardwired to 2 phenotypes) to an
# arbitrary set of continuous phenotypes. Produces:
#   1. Faceted lollipop  - top pathways per phenotype (competitive P, with the
#      permutation floor drawn explicitly; colour = PRS R2, size = N SNPs).
#   2. Pathway x phenotype heatmap - -log10(competitive P) over the union of
#      each phenotype's top pathways, to see shared vs biomarker-specific hits.
#   3. Ranked CSV table (supplement).
#
# Deps: readr, dplyr, tidyr, stringr, ggplot2 (no optparse / ggrepel).
#
# Usage (on Minerva, from the repo root, with a tidyverse-capable R):
#   Rscript scripts/plot_prset_biomarkers.R \
#     [--resultsdir <dir>] [--outdir figures/biomarkers] \
#     [--perm 10000] [--top 10]
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
outdir  <- getarg("--outdir", "figures/biomarkers")
perm    <- as.integer(getarg("--perm", "10000"))
top     <- as.integer(getarg("--top", "10"))
godict  <- getarg("--godict",
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/gene_ontology/qced_data/aux_files/GO.dict")
mgidict <- getarg("--mgidict",
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/mgi/qced_data/MGI.dict")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
floor_p    <- 1 / (perm + 1)
floor_logp <- -log10(floor_p)

# ---- phenotype panel ------------------------------------------------------
# Each row -> one PRSet .summary file (APOE excluded, full EUR sample).
# `file` is the summary basename; biomarkers follow prset_nothreshold_eur_<key>,
# while AD diagnosis / age-of-onset use their own out-prefixes.
pheno_tbl <- data.frame(
  label = c("AD diagnosis", "Age of onset", "GFAP", "NEFL",
            "Hippocampal volume", "WMH volume"),
  group = c("AD phenotype", "AD phenotype", "Plasma protein", "Plasma protein",
            "MRI", "MRI"),
  file  = c("ad_case.control_prset_nothreshold_eur.summary",
            "ad_ageofonset_prset_nothreshold_eur.summary",
            "prset_nothreshold_eur_GFAP.summary",
            "prset_nothreshold_eur_NEFL.summary",
            "prset_nothreshold_eur_Hipp_Mean_Vol.summary",
            "prset_nothreshold_eur_WMH_Vol.summary"),
  stringsAsFactors = FALSE
)
pheno_tbl$path <- file.path(resultsdir, pheno_tbl$file)
# display order: AD phenotypes first, then proteins, then MRI, in table order
pheno_levels <- pheno_tbl$label

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
read_prset <- function(path, label, group) {
  if (!file.exists(path)) { message("MISSING: ", path); return(NULL) }
  df <- read_table(path, show_col_types = FALSE)
  r2_col <- if ("PRS.R2" %in% names(df)) "PRS.R2"
            else if ("R2" %in% names(df)) "R2" else NA_character_
  if (is.na(r2_col)) stop("No R2 column in ", path)
  if (!"Competitive.P" %in% names(df))
    stop("No Competitive.P column in ", path, " (was --set-perm used?)")
  df %>%
    rename(R2 = all_of(r2_col)) %>%
    mutate(Phenotype = label, Group = group) %>%
    filter(!Set %in% c("Base", "Background"))
}

dat <- bind_rows(Map(read_prset, pheno_tbl$path, pheno_tbl$label, pheno_tbl$group))
if (is.null(dat) || nrow(dat) == 0) stop("No PRSet summaries could be read.")
dat$Phenotype <- factor(dat$Phenotype, levels = pheno_levels)

# ---- per-phenotype ranking ------------------------------------------------
ranked <- dat %>%
  mutate(at_floor   = Competitive.P <= floor_p + 1e-12,
         logP       = -log10(pmax(Competitive.P, floor_p)),
         R2_per_SNP = R2 / Num_SNP) %>%
  group_by(Phenotype) %>%
  arrange(Competitive.P, desc(R2_per_SNP), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  slice_head(n = top) %>%
  ungroup()

# ---- Figure 1: faceted lollipop -------------------------------------------
r1 <- ranked %>% mutate(Set_ord = paste(Set, Phenotype, sep = "___"))
lvls <- r1 %>% arrange(Phenotype, desc(rank)) %>% pull(Set_ord)
r1$Set_ord <- factor(r1$Set_ord, levels = unique(lvls))
strip_ord <- function(x) prettify_set(str_remove(x, "___.*$"))

p1 <- ggplot(r1, aes(logP, Set_ord)) +
  geom_vline(xintercept = floor_logp, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = logP, yend = Set_ord),
               colour = "grey80", linewidth = 0.5) +
  geom_point(aes(colour = R2, size = Num_SNP)) +
  scale_y_discrete(labels = strip_ord) +
  scale_colour_viridis_c(option = "C", name = expression(PRS~R^2)) +
  scale_size_continuous(name = "N SNPs", range = c(2, 7)) +
  facet_wrap(~ Phenotype, scales = "free_y", ncol = 2) +
  labs(x = expression(-log[10]~"(competitive P)"), y = NULL,
       title = sprintf("Top %d pathways per phenotype (competitive P)", top),
       subtitle = "Dashed line = permutation floor (1/(perm+1))") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(outdir, "prset_biomarkers_lollipop.pdf"), p1,
       width = 13, height = 8)

# ---- Figure 2: pathway x phenotype heatmap --------------------------------
top_sets <- unique(ranked$Set)
hm <- dat %>%
  filter(Set %in% top_sets) %>%
  mutate(logP = -log10(pmax(Competitive.P, floor_p)),
         at_floor = Competitive.P <= floor_p + 1e-12) %>%
  select(Set, Phenotype, logP, at_floor)

# row order: most broadly / strongly significant pathways at top
set_order <- hm %>% group_by(Set) %>%
  summarise(tot = sum(logP, na.rm = TRUE), .groups = "drop") %>%
  arrange(tot) %>% pull(Set)
hm$Set <- factor(hm$Set, levels = set_order)
hm$cell <- ifelse(hm$at_floor, sprintf("%.1f*", hm$logP), sprintf("%.1f", hm$logP))

p2 <- ggplot(hm, aes(Phenotype, Set, fill = logP)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = cell), size = 2.5, colour = "grey15") +
  scale_fill_viridis_c(option = "D", name = expression(-log[10]~"(comp. P)")) +
  scale_y_discrete(labels = function(x) prettify_set(x)) +
  labs(x = NULL, y = NULL,
       title = "Pathway x phenotype competitive-P landscape",
       subtitle = "Union of each phenotype's top pathways; * = tied at permutation floor") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())

n_sets  <- length(set_order)
n_pheno <- length(pheno_levels)
ggsave(file.path(outdir, "prset_biomarkers_heatmap.pdf"), p2,
       width = max(8, 1.4 * n_pheno + 2),
       height = max(4, 0.28 * n_sets + 2), limitsize = FALSE)

# ---- ranked table ---------------------------------------------------------
ranked %>%
  transmute(Phenotype, Group, rank, Set, Set_pretty = prettify_set(Set),
            Competitive.P, at_floor, R2, R2_per_SNP, Coefficient, Num_SNP) %>%
  arrange(Phenotype, rank) %>%
  write_csv(file.path(outdir, "prset_biomarkers_top_pathways.csv"))

message("Done. Wrote 2 figures + prset_biomarkers_top_pathways.csv to ", outdir)
n_floor <- sum(ranked$at_floor)
if (n_floor > 0)
  message(sprintf("Note: %d of the top pathways are tied at the permutation floor; ", n_floor),
          "re-run the top tier with a larger --set-perm to break ties.")
