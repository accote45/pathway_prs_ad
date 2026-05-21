library(data.table)
library(tidyverse)
library(readxl)

excel_path <- "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/coexpression/snRNAseq_DLPFC_modules/raw_data/40478_2025_2143_MOESM3_ESM.xlsx"

# Discover which sheets contain module assignment data
all_sheets <- excel_sheets(excel_path)

expected_cols <- c("ensembl", "gene_name", "module_assignment", "module_number")

module_data <- map_dfr(all_sheets, function(sheet) {
  df <- read_excel(excel_path, sheet = sheet)
  if (all(expected_cols %in% tolower(colnames(df)))) {
    colnames(df) <- tolower(colnames(df))
    df %>% select(all_of(expected_cols))
  } else {
    NULL
  }
})

message("Read ", n_distinct(module_data$module_assignment), " modules across ",
        n_distinct(sub("_m\\d+$", "", module_data$module_assignment)), " cell types")

# Summarise module sizes and filter
module_summary <- module_data %>%
  group_by(module_assignment) %>%
  summarise(genes = list(ensembl), n_genes = n(), .groups = "drop") %>%
  arrange(desc(n_genes))

module_summary_filtered <- module_summary %>%
  filter(n_genes >= 10 & n_genes <= 2000)

message("Retained ", nrow(module_summary_filtered), " of ", nrow(module_summary),
        " modules after size filtering (", nrow(module_summary) - nrow(module_summary_filtered), " removed)")
message("Filtered module size distribution:")
print(summary(module_summary_filtered$n_genes))

# Build GMT lines: name \t description \t gene1 \t gene2 \t ...
gmt_lines <- module_summary_filtered %>%
  rowwise() %>%
  mutate(
    line = paste(
      c(module_assignment, "PLACEHOLDER", unlist(genes)), collapse = "\t"
    )
  ) %>%
  pull(line)

out_path <- "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/coexpression/snRNAseq_DLPFC_modules/snRNAseq_DLPFC_modules.gmt"
writeLines(gmt_lines, out_path)
message("Written ", length(gmt_lines), " gene sets to: ", out_path)








