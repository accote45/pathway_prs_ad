

library(tidyverse)
library(data.table)
library(readxl)

excel_path <- "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/coexpression/phg_protein_modules/raw_data/" 

# Read the third tab, which contains: Protein, Symbol, is.hub, module
df <- read_excel(excel_path, sheet = 3)
colnames(df) <- tolower(colnames(df))

# Summarise module sizes and filter
module_summary <- df %>%
  group_by(module) %>%
  summarise(genes = list(symbol), n_genes = n(), .groups = "drop") %>%
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
  mutate(line = paste(c(module, "na", unlist(genes)), collapse = "\t")) %>%
  pull(line)

out_path <- "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/coexpression/phg_protein_modules/qced_data/phg_proteins.gmt"  
writeLines(gmt_lines, out_path)
message("Written ", length(gmt_lines), " gene sets to: ", out_path)












