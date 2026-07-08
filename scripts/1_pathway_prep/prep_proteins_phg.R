

library(tidyverse)
library(data.table)
library(readxl)

excel_path <- paste0(
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/",
  "coexpression/phg_protein_modules/raw_data/mmc5.xlsx"
)
ensembl_path <- paste0(
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/",
  "ensembl/qc/Ensembl.regions"
)

# Read the third tab, which contains: Protein, Symbol, is.hub, module
df <- read_excel(excel_path, sheet = 3)
colnames(df) <- tolower(colnames(df))

ensembl <- fread(ensembl_path)

get_pathway_string <- function(name, genes) {
  gene_str <- paste(genes, sep = "\t", collapse = "\t")
  with_desc <- paste("PLACEHOLDER", gene_str, sep = "\t", collapse = "\t")
  paste(name, with_desc, sep = "\t", collapse = "\t")
}

# Summarise module sizes and filter
module_summary <- df %>%
  group_by(module) %>%
  summarise(genes = list(symbol), n_genes = n(), .groups = "drop") %>%
  arrange(desc(n_genes))

module_summary_filtered <- module_summary %>%
  filter(n_genes >= 10 & n_genes <= 2000)

message("Retained ", nrow(module_summary_filtered), " of ", nrow(module_summary),
        " modules after size filtering (",
        nrow(module_summary) - nrow(module_summary_filtered), " removed)")
message("Filtered module size distribution:")
print(summary(module_summary_filtered$n_genes))

out_path <- paste0(
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/",
  "coexpression/phg_protein_modules/qced_data/phg_proteins.gmt"
)
file_conn <- file(out_path, open = "wt")

n_written <- 0L
for (i in seq_len(nrow(module_summary_filtered))) {
  row          <- module_summary_filtered[i, ]
  gene_symbols <- unique(trimws(unlist(row$genes)))
  gene_symbols <- gene_symbols[nchar(gene_symbols) > 0]

  ensembl_genes <- unique(ensembl[Name %in% gene_symbols, ID])

  if (length(ensembl_genes) >= 10 && length(ensembl_genes) <= 2000) {
    writeLines(get_pathway_string(row$module, ensembl_genes), file_conn)
    n_written <- n_written + 1L
  }
}

close(file_conn)
message("Written ", n_written, " gene sets to: ", out_path)
