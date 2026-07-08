library(data.table)
library(tidyverse)
library(readxl)

setwd('/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/synGO/')

path <- "./raw_data/ontologies.xlsx"

syngo_pathways <- read_excel(path)

ensembl <- fread("/sc/arion/projects/paul_oreilly/data/Functional_Genomics/ensembl/qc/Ensembl.regions")

cur.gmt <- "./qced_data/SynGO.gmt"
fileConn <- file(cur.gmt, open = "wt")

get_pathway_string <- function(name, genes) {
  paste(genes, sep = "\t", collapse = "\t") %>%
    paste("PLACEHOLDER", ., sep = "\t", collapse = "\t") %>%
    paste(name, ., sep = "\t", collapse = "\t") %>%
    return
}

for (i in seq_len(nrow(syngo_pathways))) {
  pathway_name <- syngo_pathways$name[i]

  gene_symbols <- syngo_pathways$hgnc_symbol[i] %>%
    strsplit(",") %>%
    unlist() %>%
    trimws() %>%
    unique()

  ensembl_genes <- ensembl[Name %in% gene_symbols, ID] %>% unique()

  if (length(ensembl_genes) >= 10 & length(ensembl_genes) <= 2000) {
    get_pathway_string(pathway_name, ensembl_genes) %>%
      writeLines(., fileConn)
  }
}

close(fileConn)