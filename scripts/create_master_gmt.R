## create master GMT file


library(tidyverse)
library(data.table)

# Read celltype specificity GMT (no description/placeholder field)
celltype_lines <- readLines("/sc/arion/projects/psychgen/cotea02_prset/judit_revisions/software/1kg_test/GeneExpressionLandscape/intermediate_files/celltype_specificity.gmt")

# Add PLACEHOLDER as second field to match standard GMT format: name\tdescription\tgenes...
celltype_lines_fixed <- sapply(celltype_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  paste(c(parts[1], "PLACEHOLDER", parts[-1]), collapse = "\t")
}, USE.NAMES = FALSE)

# Read combined pathways GMT (already has PLACEHOLDER as second field)
combined_lines <- readLines("/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/scripts/output/pathways_retained.gmt")

# Combine and write master GMT
master_lines <- c(combined_lines, celltype_lines_fixed)

# Replace spaces or tabs in pathway names (first field) with underscores
master_lines <- sapply(master_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  parts[1] <- gsub("[ \t]", "_", parts[1])
  paste(parts, collapse = "\t")
}, USE.NAMES = FALSE)

out_path <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/scripts/output/master.gmt"
writeLines(master_lines, out_path)

cat("Master GMT created:", length(master_lines), "pathways\n")
cat("  - Combined pathways:", length(combined_lines), "\n")
cat("  - Cell type specificity:", length(celltype_lines), "\n")

# Create APOE-excluded version
apoe_gene <- "ENSG00000130203"

master_lines_no_apoe <- sapply(master_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  parts <- parts[parts != apoe_gene]
  paste(parts, collapse = "\t")
}, USE.NAMES = FALSE)

out_path_no_apoe <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/scripts/output/master_no_apoe.gmt"
writeLines(master_lines_no_apoe, out_path_no_apoe)

cat("APOE-excluded GMT created:", out_path_no_apoe, "\n")

# Extract all unique genes from master GMT for PRSet background (APOE excluded)
all_genes <- unique(unlist(lapply(master_lines_no_apoe, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  parts[-(1:2)]  # drop name and description fields
})))

background_path <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/scripts/output/background_genes.txt"
writeLines(all_genes, background_path)

cat("Background genes file created:", length(all_genes), "unique genes\n")




