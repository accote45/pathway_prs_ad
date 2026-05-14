library(tidyverse)
library(data.table)

# --- Load summary files ---
master <- fread("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/prset/birewire/msigdbgenes/bmi/bmi_set.birewire.summary", sep = "\t")
eur    <- fread("/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/bmi/bmi_prset_nothreshold_eur.summary", sep = "\t")
afr    <- fread("/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/bmi/bmi_prset_nothreshold_afr.summary", sep = "\t")

# --- Reference universe: pathways positive in master ---
master_pos <- master[Set != "Base" & !is.na(Competitive.P)] %>%
  select(Set, Competitive.P_Master = Competitive.P, PRS.R2_Master = PRS.R2,
         Coefficient_Master = Coefficient, Num_SNP_Master = Num_SNP)

# --- EUR and AFR: no positivity filter, just require non-NA Competitive.P ---
eur_data <- eur[Set != "Base" & !is.na(Competitive.P)] %>%
  select(Set, Competitive.P_EUR = Competitive.P, PRS.R2_EUR = PRS.R2,
         Coefficient_EUR = Coefficient, Num_SNP_EUR = Num_SNP)

afr_data <- afr[Set != "Base" & !is.na(Competitive.P)] %>%
  select(Set, Competitive.P_AFR = Competitive.P, PRS.R2_AFR = PRS.R2,
         Coefficient_AFR = Coefficient, Num_SNP_AFR = Num_SNP)

# --- Join EUR and AFR to master-positive pathways ---
all_pos <- master_pos %>%
  inner_join(eur_data, by = "Set") %>%
  inner_join(afr_data, by = "Set")

# --- Rank all three within the shared set ---
all_pos <- all_pos %>%
  mutate(
    Ratio_Master = PRS.R2_Master / Num_SNP_Master,
    Ratio_EUR    = PRS.R2_EUR    / Num_SNP_EUR,
    Ratio_AFR    = PRS.R2_AFR    / Num_SNP_AFR
  ) %>%
  arrange(Competitive.P_Master, desc(Ratio_Master)) %>%
  mutate(Rank_Master = row_number()) %>%
  arrange(Competitive.P_EUR, desc(Ratio_EUR)) %>%
  mutate(Rank_EUR = row_number()) %>%
  arrange(Competitive.P_AFR, desc(Ratio_AFR)) %>%
  mutate(Rank_AFR = row_number()) %>%
  arrange(Rank_Master)

cat("Pathways positive in master with EUR and AFR data:", nrow(all_pos), "\n\n")

# --- Combined table ---
cat("=== Ranked by Master Competitive.P ===\n")
print(all_pos %>%
  select(Set, Competitive.P_Master, Rank_Master,
         Competitive.P_EUR, Rank_EUR,
         Competitive.P_AFR, Rank_AFR,
         Coefficient_Master, Coefficient_EUR, Coefficient_AFR,
         Num_SNP_Master, Num_SNP_EUR, Num_SNP_AFR))

# --- Save ---
fwrite(all_pos, "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/bmi/crossancestry_pathway_comparison.tsv", sep = "\t")

# --- Kendall rank correlation ---
run_kendall <- function(df, rank_x, rank_y, label) {
  result <- cor.test(df[[rank_x]], df[[rank_y]], method = "kendall")
  cat(sprintf("%s (n=%d): tau = %.4f, z = %.4f, p = %.4g\n",
              label, nrow(df), result$estimate, result$statistic, result$p.value))
}

cat("\n=== Kendall Rank Correlation (Master vs EUR) ===\n")
run_kendall(all_pos,                                  "Rank_Master", "Rank_EUR", "All shared pathways")
run_kendall(all_pos %>% filter(Rank_Master <= 20),   "Rank_Master", "Rank_EUR", "Top 20 Master")
run_kendall(all_pos %>% filter(Rank_Master <= 50),   "Rank_Master", "Rank_EUR", "Top 50 Master")
run_kendall(all_pos %>% filter(Rank_Master <= 100),  "Rank_Master", "Rank_EUR", "Top 100 Master")
run_kendall(all_pos %>% filter(Rank_Master <= 500),  "Rank_Master", "Rank_EUR", "Top 500 Master")
run_kendall(all_pos %>% filter(Rank_Master <= 1000), "Rank_Master", "Rank_EUR", "Top 1000 Master")

cat("\n=== Kendall Rank Correlation (Master vs AFR) ===\n")
run_kendall(all_pos,                                  "Rank_Master", "Rank_AFR", "All shared pathways")
run_kendall(all_pos %>% filter(Rank_Master <= 20),   "Rank_Master", "Rank_AFR", "Top 20 Master")
run_kendall(all_pos %>% filter(Rank_Master <= 50),   "Rank_Master", "Rank_AFR", "Top 50 Master")
run_kendall(all_pos %>% filter(Rank_Master <= 100),  "Rank_Master", "Rank_AFR", "Top 100 Master")
run_kendall(all_pos %>% filter(Rank_Master <= 500),  "Rank_Master", "Rank_AFR", "Top 500 Master")
run_kendall(all_pos %>% filter(Rank_Master <= 1000), "Rank_Master", "Rank_AFR", "Top 1000 Master")

# --- Scatterplots: Top 50 and Top 20 Master, vs EUR and AFR ---
top50 <- all_pos %>% filter(Rank_Master <= 50)
top20 <- all_pos %>% filter(Rank_Master <= 20)

make_scatter <- function(df, y_rank, y_p, y_label, title_suffix, outfile) {
  k <- cor.test(df$Rank_Master, df[[y_rank]], method = "kendall")
  ggplot(df, aes(x = Rank_Master, y = .data[[y_rank]], label = Set)) +
    geom_point(aes(color = Competitive.P_Master), size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
    geom_text(size = 2.5, vjust = -0.8, hjust = 0.5, check_overlap = TRUE) +
    scale_color_viridis_c(name = "Master\nCompetitive.P", direction = -1) +
    labs(
      title = paste0("Top ", nrow(df), " Master Pathways: Master vs ", title_suffix),
      subtitle = sprintf("Kendall tau = %.4f, p = %.4g (n=%d)",
                         k$estimate, k$p.value, nrow(df)),
      x = "Master Rank",
      y = paste0(y_label, " Rank (within shared pathways)")
    ) +
    theme_bw()
  ggsave(outfile, width = 8, height = 7)
}

make_scatter(top50, "Rank_EUR", "Competitive.P_EUR", "EUR", "EUR Rank", "top50_master_eur_rank_scatter.pdf")
make_scatter(top20, "Rank_EUR", "Competitive.P_EUR", "EUR", "EUR Rank", "top20_master_eur_rank_scatter.pdf")
make_scatter(top50, "Rank_AFR", "Competitive.P_AFR", "AFR", "AFR Rank", "top50_master_afr_rank_scatter.pdf")
make_scatter(top20, "Rank_AFR", "Competitive.P_AFR", "AFR", "AFR Rank", "top20_master_afr_rank_scatter.pdf")