#!/usr/bin/env Rscript
# Combine PRSet .summary files into one multi-tab Excel workbook.
#
# Per input file:
#   * drop the genome-wide Base / Background sets (keep pathways only)
#   * add a Pathway_Name column: GO:####/MP:#### IDs are resolved to their
#     term names via the GO/MGI dicts; already-named sets (e.g. MSigDB) are
#     tidied but keep their name
#   * sort by Competitive.P (ascending), then PRS.R2 / Num.SNP (descending).
#     Only PRS.R2 is ever used here — never the full-model R2.
# One worksheet per input file; tab name is derived from the file name.
#
# Usage:
#   Rscript combine_prset_summaries_xlsx.R -o prset_pathways.xlsx file1.summary file2.summary ...
#   Rscript combine_prset_summaries_xlsx.R -o out.xlsx /path/to/results/*.summary

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tools)
  library(openxlsx)
})

opt_parser <- OptionParser(option_list = list(
  make_option(c("-o", "--out"), type = "character",
              default = "prset_pathways.xlsx",
              help = "Output .xlsx path [default %default]"),
  make_option("--godict", type = "character",
              default = "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/gene_ontology/qced_data/aux_files/GO.dict",
              help = "GO ID -> term name dictionary (tab-separated: ID, Name)"),
  make_option("--mgidict", type = "character",
              default = "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/mgi/qced_data/MGI.dict",
              help = "MGI (MP) ID -> term name dictionary (tab-separated: ID, Name)")
))
opt  <- parse_args(opt_parser, positional_arguments = TRUE)
files <- opt$args
if (length(files) == 0) stop("Provide at least one .summary file.")

# ---- ID -> term name lookup (GO/MGI dicts; optional) ----------------------
# Both dicts are tab-separated with columns ID, Name; GO uses "GO:####" IDs
# and MGI uses "MP:####" IDs, so a single combined map is unambiguous.
load_dict <- function(path, label) {
  if (!file.exists(path)) {
    message(label, " dictionary not found (", path,
            "); those terms will show IDs only.")
    return(character(0))
  }
  d <- readr::read_tsv(path, show_col_types = FALSE,
                       col_types = readr::cols(.default = "c"))
  setNames(d$Name, d$ID)
}
term_map <- c(load_dict(opt$options$godict, "GO"),
              load_dict(opt$options$mgidict, "MGI"))

# Turn "MSigDB__BORCZUK_MALIGNANT_MESOTHELIOMA_UP" into
# "MSigDB: Borczuk Malignant Mesothelioma Up", "GO__GO:0043032" into
# "GO: mitochondrion inheritance", "MGI__MP:0000547" into "MGI: <name>".
prettify_set <- function(set) {
  parts <- str_split_fixed(set, "__", 2)
  db    <- parts[, 1]
  name  <- parts[, 2]
  no_sep <- name == ""                 # no "__" separator present
  db[no_sep]   <- ""
  name[no_sep] <- set[no_sep]
  # Pure ID terms ("GO:0043032" / "MP:0000547"): swap the ID for its name.
  pure_id <- str_detect(name, "^(GO|MP):[0-9]+$")
  hit     <- pure_id & name %in% names(term_map)
  name[hit] <- unname(term_map[name[hit]])
  name <- str_replace_all(name, "_", " ")
  msig <- db == "MSigDB"               # MSigDB names are ALL CAPS
  name[msig] <- str_to_title(name[msig])
  ifelse(db == "", name, paste0(db, ": ", name))
}

# Excel sheet names: <=31 chars, no : \ / ? * [ ] and must be unique.
sheet_name <- function(paths) {
  nm <- file_path_sans_ext(basename(paths))
  nm <- gsub("[:\\\\/?*\\[\\]]", "_", nm)
  nm <- substr(nm, 1, 31)
  make.unique(nm, sep = "_")
}

read_annotated <- function(path) {
  df <- readr::read_table(path, show_col_types = FALSE)

  # Only the pathway PRS.R2 is ever used for ranking — never the full-model R2.
  if (!"PRS.R2" %in% names(df)) stop("No PRS.R2 column found in ", path)
  if (!"Competitive.P" %in% names(df)) stop("No Competitive.P column in ", path)
  if (!"Num_SNP" %in% names(df)) stop("No Num_SNP column in ", path)
  if (!"Set" %in% names(df)) stop("No Set column in ", path)

  df %>%
    filter(!Set %in% c("Base", "Background")) %>%
    mutate(Pathway_Name = prettify_set(Set), .after = Set) %>%
    arrange(.data[["Competitive.P"]],
            desc(.data[["PRS.R2"]] / .data[["Num_SNP"]]))
}

wb   <- createWorkbook()
tabs <- sheet_name(files)
hdr  <- createStyle(textDecoration = "bold")

for (i in seq_along(files)) {
  d <- read_annotated(files[i])
  addWorksheet(wb, tabs[i])
  writeData(wb, tabs[i], d, headerStyle = hdr)
  freezePane(wb, tabs[i], firstRow = TRUE)
  setColWidths(wb, tabs[i], cols = seq_along(d), widths = "auto")
  message(sprintf("  %-31s %d pathways", tabs[i], nrow(d)))
}

saveWorkbook(wb, opt$options$out, overwrite = TRUE)
message("Wrote ", opt$options$out, " (", length(files), " tabs)")
