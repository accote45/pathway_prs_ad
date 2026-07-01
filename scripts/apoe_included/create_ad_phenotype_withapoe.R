# ---------------------------------------------------------------------------
# create AD case/control phenotype residuals -- APOE-INCLUDED version
#
# This is the APOE-INCLUDED counterpart of scripts/create_ad_phenotype.R.
# The ONLY analytical difference is the covariate set used to residualize AD:
# APOE_genotype is NOT regressed out here, so the APOE effect is retained in
# the phenotype (and, downstream, is allowed to contribute to the PRS via the
# APOE-containing pathway file and the removal of the chr19 --x-range).
#
#   Original (APOE regressed out): Sex + Age + Age2 + Batch + APOE_genotype + PC1-15
#   This script (APOE retained)  : Sex + Age + Age2 + Batch + PC1-15
#
# Output: ad_phenotype_residuals_withapoe.txt   (residual column: AD_resid)
# The residual column keeps the name AD_resid so the same --pheno-col works.
# ---------------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(dplyr)

data_dir <- '/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data'

# Read in all files
ad_dementia   <- fread('/sc/arion/projects/paul_oreilly/lab/shared/pheno/AD_Dementia_28March2024.txt')
apoe_geno     <- fread('/sc/arion/projects/paul_oreilly/lab/shared/UKB_allSample_APOE_genotype.txt')
alz_case_ctrl <- fread('/sc/arion/projects/paul_oreilly/lab/shared/pheno/clive/alz_case_control.csv')
covs          <- fread('/sc/arion/projects/paul_oreilly/data/ukb/phenotype/ukb18177.covar')
age_sex_covar <- fread('/sc/arion/projects/paul_oreilly/lab/shared/pheno/ageSex.covar')

dat <- ad_dementia %>%
  inner_join(apoe_geno, by = c("sample_id" = "IID")) %>%
  left_join(alz_case_ctrl[, .(sample_id, dob, death_date, alz_diag_date)], by = "sample_id") %>%
  inner_join(covs, by = c("sample_id" = "IID")) %>%
  inner_join(age_sex_covar, by = c("sample_id" = "IID"))

eur_id_file <- file.path(data_dir, 'eur_sample_ids_80pc.txt')
eur_ids <- fread(eur_id_file, header = FALSE)$V1
dat <- dat[dat$sample_id %in% eur_ids, ]
# e1 alleles are still dropped (rare / QC), matching the original pipeline.
dat <- dat[!grepl("e1", dat$APOE_genotype, fixed = TRUE), ]

dat$Age2 <- dat$Age^2
dat$APOE_genotype[dat$APOE_genotype == ""] <- NA
dat$APOE_genotype <- as.factor(dat$APOE_genotype)

table(dat$AD)

# ---------------------------------------------------------------------------
# Adjust AD diagnosis for covariates, obtain residuals
# APOE-INCLUDED: APOE_genotype is deliberately excluded from the covariate set.
# ---------------------------------------------------------------------------

covariates <- c("Sex", "Age", "Age2", "Batch",
                paste0("PC", 1:15))

covariate_formula <- paste(covariates, collapse = " + ")

# For binary outcomes, use logistic regression with Pearson residuals
get_logistic_residuals <- function(outcome, data) {
  model_cols <- c(outcome, covariates)
  complete_rows <- complete.cases(data[, model_cols, with = FALSE])
  complete_data <- data[complete_rows, ]

  formula_str <- paste(outcome, "~", covariate_formula)
  model <- glm(as.formula(formula_str), data = complete_data, family = binomial)

  residuals_vector <- rep(NA, nrow(data))
  residuals_vector[complete_rows] <- residuals(model, type = "pearson")

  return(residuals_vector)
}

binary_outcomes <- c("AD")

final_pheno <- dat
for (outcome in binary_outcomes) {
  if (outcome %in% colnames(final_pheno)) {
    cat("Adjusting", outcome, "for covariates (APOE NOT regressed out)...\n")
    final_pheno[[paste0(outcome, "_resid")]] <- get_logistic_residuals(outcome, final_pheno)
  }
}

# Save the adjusted phenotype file (APOE-included residual)
out_path <- file.path(data_dir, 'ad_phenotype_residuals_withapoe.txt')
write.table(final_pheno, file = out_path, sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("\nWrote:", out_path, "\n")
resid_cols <- colnames(final_pheno)[grepl("_resid$", colnames(final_pheno))]
cat("Adjusted outcomes available:\n", paste(resid_cols, collapse = "\n"), "\n")
