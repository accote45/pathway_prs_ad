# create master phenotype file
# read in case/control status + age of onset
# read in covariates (age, sex, PCs)
# get phenotype residuals


# ancestry assignments: /sc/arion/projects/paul_oreilly/data/Biobanks/CrossBiobank_qc/GeneticDerived_Ancestry/admixture/data/UKB/ukb_ancestry_assignment.txt
  ## will use ADMIXTURE proportion threshold of 80% to define AFR ancestry individual
  n=7,170

# create ancestry assignment files
awk 'NR>1 && $4=="AFR" {print $2}' /sc/arion/projects/paul_oreilly/data/Biobanks/CrossBiobank_qc/GeneticDerived_Ancestry/admixture/data/UKB/ukb_ancestry_assignment.txt > /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/afr_sample_ids_80pc.txt

awk 'NR>1 && $4=="EUR" {print $2}' /sc/arion/projects/paul_oreilly/data/Biobanks/CrossBiobank_qc/GeneticDerived_Ancestry/admixture/data/UKB/ukb_ancestry_assignment.txt > /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/eur_sample_ids_80pc.txt



# file purgatory:
    # don't need, same onset date as other file: /sc/arion/projects/paul_oreilly/lab/shared/pheno/clive/Alzheimer_disease_case.txt

    /sc/arion/projects/paul_oreilly/lab/shared/pheno/AD_Dementia_28March2024.txt
    /sc/arion/projects/paul_oreilly/lab/shared/UKB_allSample_APOE_genotype.txt
    /sc/arion/projects/paul_oreilly/lab/shared/pheno/ageSex.covar
    /sc/arion/projects/paul_oreilly/lab/shared/pheno/clive/alz_case_control.csv
    /sc/arion/projects/paul_oreilly/lab/shared/pheno/clive/Alzheimer_disease_case.txt



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
dat <- dat[!grepl("e1", dat$APOE_genotype, fixed = TRUE), ]

dat$Age2 <- dat$Age^2
dat$APOE_genotype[dat$APOE_genotype == ""] <- NA
dat$APOE_genotype <- as.factor(dat$APOE_genotype)

table(dat$AD)

# Basic demographics by AD status
dat %>%
  group_by(AD) %>%
  summarise(
    n = n(),

    # Age
    age_mean = mean(Age, na.rm = TRUE),
    age_sd   = sd(Age, na.rm = TRUE),

    # Sex (assuming 0 = male, 1 = female)
    n_female = sum(Sex == 1, na.rm = TRUE),
    pct_female = mean(Sex == 1, na.rm = TRUE) * 100,

    # APOE e4
    n_e4_carrier  = sum(e4_copy >= 1, na.rm = TRUE),
    pct_e4_carrier = mean(e4_copy >= 1, na.rm = TRUE) * 100,
    n_e4_homozygous = sum(e4_copy == 2, na.rm = TRUE),

    # APOE e2
    n_e2_carrier  = sum(e2_copy >= 1, na.rm = TRUE),
    pct_e2_carrier = mean(e2_copy >= 1, na.rm = TRUE) * 100,

    .groups = "drop"
  )

# APOE genotype breakdown
table(dat$APOE_genotype, dat$AD)

# Chi-square tests for categorical variables
chisq.test(table(dat$Sex, dat$AD))
chisq.test(table(dat$e4_copy, dat$AD))

# T-test for age
t.test(Age ~ AD, data = dat)


# Basic demographics by case/control status
get_demographics <- function(data, outcome_col) {
  data %>%
    group_by(.data[[outcome_col]]) %>%
    summarise(
      outcome = outcome_col,
      n = n(),

      # Age
      age_mean = mean(Age, na.rm = TRUE),
      age_sd   = sd(Age, na.rm = TRUE),

      # Sex (0 = male, 1 = female)
      n_female   = sum(Sex == 1, na.rm = TRUE),
      pct_female = mean(Sex == 1, na.rm = TRUE) * 100,

      # APOE e4
      n_e4_carrier    = sum(e4_copy >= 1, na.rm = TRUE),
      pct_e4_carrier  = mean(e4_copy >= 1, na.rm = TRUE) * 100,
      n_e4_homozygous = sum(e4_copy == 2, na.rm = TRUE),

      # APOE e2
      n_e2_carrier   = sum(e2_copy >= 1, na.rm = TRUE),
      pct_e2_carrier = mean(e2_copy >= 1, na.rm = TRUE) * 100,

      .groups = "drop"
    )
}

outcomes <- c("AD", "VD", "Dementia")
demo_table <- bind_rows(lapply(outcomes, get_demographics, data = dat))
print(demo_table)

# APOE genotype breakdown by outcome
for (outcome in outcomes) {
  cat("\n--- APOE genotype x", outcome, "---\n")
  print(table(dat$APOE_genotype, dat[[outcome]]))
}



# ---------------------------------------------------------------------------
# Adjust AD diagnosis for covariates, obtain residuals
# ---------------------------------------------------------------------------

# Define covariates (includes APOE genotype as factor)
covariates <- c("Sex", "Age", "Age2", "Batch",
                "APOE_genotype",
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

# Apply logistic regression adjustment to AD
binary_outcomes <- c("AD")

final_pheno <- dat
for (outcome in binary_outcomes) {
  if (outcome %in% colnames(final_pheno)) {
    cat("Adjusting", outcome, "for covariates using logistic regression...\n")
    final_pheno[[paste0(outcome, "_resid")]] <- get_logistic_residuals(outcome, final_pheno)
  }
}

# Save the adjusted phenotype file
out_path <- file.path(data_dir, 'ad_phenotype_residuals.txt')
write.table(final_pheno, file = out_path, sep = "\t",
row.names = FALSE, quote = FALSE)

# Print summary of adjustments
cat("\nAdjusted outcomes available:\n")
resid_cols <- colnames(final_pheno)[grepl("_resid$", colnames(final_pheno))]
cat(paste(resid_cols, collapse = "\n"))
