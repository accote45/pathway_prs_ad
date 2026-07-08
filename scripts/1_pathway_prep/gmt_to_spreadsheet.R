## Convert a GMT file into a tidy, human-readable spreadsheet.
## Output columns: Pathway | Source | Description | Num_genes | Genes (comma-separated)
##
## GMT is a machine format (ragged, unlabeled, tab-delimited). This produces a
## CSV a non-specialist can open in Excel, sort, and search.
##
## Descriptions are resolved per source (see resolve_description() below):
##   - GO terms  -> looked up from the GO.db Bioconductor package (offline)
##   - anything that already has a real description in the GMT -> kept as-is
##   - everything else -> the pathway name itself, cleaned up (it's already readable)

library(tidyverse)

gmt_path <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/pathways/master.gmt"
out_path <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/pathways/master_pathways.csv"

# Each line: name \t description \t gene1 \t gene2 \t ...  (variable # of genes)
lines <- readLines(gmt_path)

df <- map_dfr(lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  genes <- parts[-(1:2)]                 # drop name + description
  genes <- genes[genes != "" & !is.na(genes)]
  tibble(
    Pathway  = parts[1],
    RawDesc  = if (length(parts) >= 2) parts[2] else NA_character_,
    Num_genes = length(genes),
    Genes    = paste(genes, collapse = ", ")
  )
})

# Source is the prefix before the first "__" (see prune_pathways.R). Celltype
# sets added afterwards have no "__" and no prefix; label them "Celltype".
df <- df %>%
  mutate(Source = if_else(str_detect(Pathway, "__"),
                          sub("__.*", "", Pathway),
                          "Celltype"))

## --- GO term lookup ---------------------------------------------------------
## Map GO IDs (e.g. "GO:0016462") to their term names via GO.db. If the package
## isn't installed, we warn and fall back to the cleaned-up name for GO sets too.
go_pathways <- df %>% filter(Source == "GO")
go_terms <- tibble(GOID = character(), TERM = character())

if (nrow(go_pathways) > 0) {
  if (requireNamespace("GO.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    go_ids <- unique(sub("^GO__", "", go_pathways$Pathway))   # -> "GO:0016462"
    go_terms <- AnnotationDbi::select(GO.db::GO.db,
                                      keys = go_ids,
                                      columns = "TERM",
                                      keytype = "GOID") %>%
      as_tibble() %>%
      distinct(GOID, .keep_all = TRUE)
    n_missing <- sum(is.na(go_terms$TERM))
    cat("GO.db: resolved", sum(!is.na(go_terms$TERM)), "of",
        length(go_ids), "GO terms",
        if (n_missing > 0) paste0("(", n_missing, " unmatched/obsolete)") else "",
        "\n")
  } else {
    cat("NOTE: GO.db / AnnotationDbi not installed -- GO descriptions will fall",
        "back to the raw ID.\n",
        "      Install once with:",
        "BiocManager::install('GO.db')\n")
  }
}
go_lookup <- setNames(go_terms$TERM, go_terms$GOID)

# Clean a pathway name into a readable label: strip "SOURCE__", underscores -> spaces
clean_name <- function(pathway) {
  nm <- sub("^[^_]+__", "", pathway)   # drop source prefix if present
  str_squish(gsub("_", " ", nm))
}

resolve_description <- function(pathway, source, raw_desc) {
  # 1) GO terms: use the looked-up term name when available
  if (source == "GO") {
    goid <- sub("^GO__", "", pathway)
    if (!is.na(go_lookup[goid]) && nzchar(go_lookup[goid])) return(unname(go_lookup[goid]))
  }
  # 2) A real description already in the GMT (e.g. ClinPGx/PharmGKB names)
  if (!is.na(raw_desc) && nzchar(raw_desc) && raw_desc != "PLACEHOLDER") return(raw_desc)
  # 3) Fall back to the (already readable) pathway name
  clean_name(pathway)
}

df <- df %>%
  mutate(Description = pmap_chr(list(Pathway, Source, RawDesc), resolve_description)) %>%
  arrange(Source, Pathway) %>%
  select(Pathway, Source, Description, Num_genes, Genes)

# write.csv quotes the Genes field automatically (it contains commas), so the
# CSV stays valid even with thousands of genes in one cell.
write.csv(df, out_path, row.names = FALSE)

cat("Wrote", nrow(df), "pathways to", out_path, "\n")
cat("Gene counts: min", min(df$Num_genes),
    "median", median(df$Num_genes),
    "max", max(df$Num_genes), "\n")
cat("Pathways per source:\n")
print(count(df, Source, name = "n_pathways"))

# Heads-up: Excel truncates any single cell over 32,767 characters on display
# (the CSV itself is fine). Flag pathways whose gene list may hit that limit.
too_long <- df %>% filter(nchar(Genes) > 32767)
if (nrow(too_long) > 0) {
  cat("\nWARNING:", nrow(too_long),
      "pathway(s) have a gene list >32,767 chars and will be truncated in Excel:\n")
  cat(paste0("  - ", too_long$Pathway, " (", too_long$Num_genes, " genes)"),
      sep = "\n")
}
