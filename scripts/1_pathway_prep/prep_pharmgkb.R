library(dplyr)
library(stringr)
library(data.table)

ensembl <- fread("/sc/arion/projects/paul_oreilly/data/Functional_Genomics/ensembl/qc/Ensembl.regions")

parse_pharmgkb_pathway <- function(file) {
  df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   quote = "", fill = TRUE)
  
  # Pathway name from filename: strip ID prefix and extension
  fname <- tools::file_path_sans_ext(basename(file))
  pathway_name <- gsub("^PA\\d+-", "", fname)
  pathway_name <- gsub("_", " ", pathway_name)
  
  # Extract genes from the "Genes" column, split comma-separated entries
  if (!"Genes" %in% colnames(df)) return(NULL)
  raw <- df$Genes
  raw <- as.character(raw)
  raw <- raw[!is.na(raw) & raw != ""]
  symbols <- unlist(strsplit(raw, ",\\s*"))
  symbols <- str_trim(symbols)
  symbols <- str_extract(symbols, "^[A-Z0-9][A-Z0-9\\-]+")
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  
  # Convert HGNC symbols to Ensembl IDs
  ensembl_ids <- unique(ensembl[Name %in% symbols, ID])
  
  # Filter by gene set size
  if (length(ensembl_ids) < 10 || length(ensembl_ids) > 2000) return(NULL)
  list(id = fname, name = pathway_name, genes = ensembl_ids)
}

pathway_files <- list.files("/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/clinpgx", 
                             pattern = "\\.tsv$", full.names = TRUE)

# Parse all pathways
pathways <- lapply(pathway_files, parse_pharmgkb_pathway)
pathways <- Filter(function(x) !is.null(x) && length(x$genes) > 0, pathways)

write_gmt <- function(pathways, outfile) {
  lines <- sapply(pathways, function(p) {
    paste(c(p$id, p$name, p$genes), collapse = "\t")
  })
  writeLines(lines, outfile)
}

out_dir <- "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/clinpgx/qced_data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_gmt(pathways, file.path(out_dir, "pharmgkb_pathways.gmt"))
