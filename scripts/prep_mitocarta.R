library(data.table)
library(tidyverse)
library(readxl)


setwd('/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/mitocarta/')

path <- "./raw_data/Human.MitoCarta3.0.xls"

sheets <- excel_sheets(path)
all_tabs <- lapply(sheets, function(s) read_excel(path, sheet = s))
names(all_tabs) <- sheets

mito_pathways <- all_tabs[[4]]

ensembl <- fread("/sc/arion/projects/paul_oreilly/data/Functional_Genomics/ensembl/qc/Ensembl.regions")

cur.gmt <- "./qced_data/MitoCarta.gmt"
fileConn <- file(cur.gmt, open = "wt")

get_pathway_string <- function(name, genes) {
  paste(genes, sep = "\t", collapse = "\t") %>%
    paste("PLACEHOLDER", ., sep = "\t", collapse = "\t") %>%
    paste(name, ., sep = "\t", collapse = "\t") %>%
    return
}

for (i in seq_len(nrow(mito_pathways))) {
  pathway_name <- mito_pathways$MitoPathway[i]
  
  gene_symbols <- mito_pathways$Genes[i] %>%
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