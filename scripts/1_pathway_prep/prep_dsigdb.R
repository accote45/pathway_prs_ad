library(data.table)
library(tidyverse)

setwd('/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/dsigdb/')

path <- "./raw_data/DSigDB_All.gmt"

ensembl <- fread("/sc/arion/projects/paul_oreilly/data/Functional_Genomics/ensembl/qc/Ensembl.regions")

cur.gmt <- "./qced_data/DSigDB.gmt"
fileConn <- file(cur.gmt, open = "wt")

get_pathway_string <- function(name, genes) {
  paste(genes, sep = "\t", collapse = "\t") %>%
    paste("PLACEHOLDER", ., sep = "\t", collapse = "\t") %>%
    paste(name, ., sep = "\t", collapse = "\t") %>%
    return
}

lines <- readLines(path)

for (line in lines) {
  fields <- strsplit(line, "\t")[[1]]
  pathway_name <- fields[1]

  # fields[2] is the URL — skip it; genes start at index 3
  gene_symbols <- fields[-(1:2)] %>%
    trimws() %>%
    unique()
  gene_symbols <- gene_symbols[nchar(gene_symbols) > 0]

  ensembl_genes <- ensembl[Name %in% gene_symbols, ID] %>% unique()

  if (length(ensembl_genes) >= 10 & length(ensembl_genes) <= 2000) {
    get_pathway_string(pathway_name, ensembl_genes) %>%
      writeLines(., fileConn)
  }
}

close(fileConn)