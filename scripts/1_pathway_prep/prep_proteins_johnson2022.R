
library(tidyverse)
library(data.table)
library(readxl)

excel_path <- paste0(
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/",
  "coexpression/johnson2022_protein_modules/raw_data/41593_2021_999_MOESM3_ESM.xlsx"
)
ensembl_path <- paste0(
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/",
  "ensembl/qc/Ensembl.regions"
)

df <- read_excel(excel_path, sheet = 5, skip = 2)
colnames(df) <- tolower(colnames(df))
print(colnames(df))

# Second column contains "SYMBOL|UniProtID" — extract the gene symbol
id_col <- colnames(df)[2]
df <- df %>%
  mutate(gene_symbol = sub("\\|.*$", "", .data[[id_col]]))

ensembl <- fread(ensembl_path)

# Update these if column names differ from what's printed above
gene_col   <- "gene_symbol"
module_col <- "net$colors"

get_pathway_string <- function(name, genes) {
  gene_str <- paste(genes, sep = "\t", collapse = "\t")
  with_desc <- paste("PLACEHOLDER", gene_str, sep = "\t", collapse = "\t")
  paste(name, with_desc, sep = "\t", collapse = "\t")
}

module_summary <- df %>%
  group_by(.data[[module_col]]) %>%
  summarise(
    genes = list(.data[[gene_col]]), n_genes = n(), .groups = "drop"
  ) %>%
  rename(module = 1) %>%
  arrange(desc(n_genes))

module_summary_filtered <- module_summary %>%
  filter(n_genes >= 10 & n_genes <= 2000)

n_removed <- nrow(module_summary) - nrow(module_summary_filtered)
n_kept  <- nrow(module_summary_filtered)
n_total <- nrow(module_summary)
message("Retained ", n_kept, " of ", n_total,
        " modules after size filtering (", n_removed, " removed)")
message("Filtered module size distribution:")
print(summary(module_summary_filtered$n_genes))

out_path <- paste0(
  "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/",
  "coexpression/johnson2022_protein_modules/qced_data/johnson2022_proteins.gmt"
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