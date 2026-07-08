# ---------------------------------------------------------------------------
# Build the APOE-INCLUDED PRSet background gene list.
#
# create_master_gmt.R already writes:
#   master.gmt          -> pathways WITH APOE (ENSG00000130203)
#   master_no_apoe.gmt  -> APOE removed
#   background_genes.txt -> background built from master_no_apoe.gmt (APOE absent)
#
# The APOE-included PRSet run uses master.gmt as --msigdb, so its --background
# must also contain APOE. This script regenerates the background from master.gmt
# so the APOE gene is present in the competitive-null background.
#
# Output: background_genes_withapoe.txt
# ---------------------------------------------------------------------------

master_path     <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/pathways/master.gmt"
background_path <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/pathways/background_genes_withapoe.txt"

master_lines <- readLines(master_path)

# Drop name (col1) and description/placeholder (col2); collect unique genes.
all_genes <- unique(unlist(lapply(master_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  parts[-(1:2)]
})))
all_genes <- all_genes[nchar(all_genes) > 0]

writeLines(all_genes, background_path)

cat("APOE-included background written:", background_path, "\n")
cat("  unique genes:", length(all_genes), "\n")
cat("  APOE (ENSG00000130203) present:",
    "ENSG00000130203" %in% all_genes, "\n")
