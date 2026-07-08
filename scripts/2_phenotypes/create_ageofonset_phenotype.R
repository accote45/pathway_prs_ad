# create master phenotype file
# read in case/control status + age of onset
# read in covariates (age, sex, PCs)
# get phenotype residuals


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

eur_id_file <- file.path(data_dir, 'samples', 'eur_sample_ids_80pc.txt')
eur_ids <- fread(eur_id_file, header = FALSE)$V1
dat <- dat[dat$sample_id %in% eur_ids, ]
dat <- dat[!grepl("e1", dat$APOE_genotype, fixed = TRUE), ]

dat$Age2 <- dat$Age^2
dat$APOE_genotype[dat$APOE_genotype == ""] <- NA
dat$APOE_genotype <- as.factor(dat$APOE_genotype)

# ---------------------------------------------------------------------------
# Calculate age of onset (years from dob to alz_diag_date)
# ---------------------------------------------------------------------------
dat$dob           <- as.Date(dat$dob)
dat$alz_diag_date <- as.Date(dat$alz_diag_date)
dat$age_of_onset  <- as.numeric(difftime(dat$alz_diag_date, dat$dob, units = "days")) / 365.25

# Filter to AD cases with a valid age of onset
cases <- dat[dat$AD == 1 & !is.na(dat$age_of_onset), ]
cat("AD cases with age of onset:", nrow(cases), "\n")
print(summary(cases$age_of_onset))

# ---------------------------------------------------------------------------
# Adjust age of onset for covariates, obtain residuals (linear regression)
# Age and Age2 excluded — age of onset is the outcome
# ---------------------------------------------------------------------------

covariates <- c("Sex", "Batch",
                "APOE_genotype",
                paste0("PC", 1:15))

covariate_formula <- paste(covariates, collapse = " + ")

get_linear_residuals <- function(outcome, data) {
  model_cols <- c(outcome, covariates)
  complete_rows <- complete.cases(data[, model_cols, with = FALSE])
  complete_data <- data[complete_rows, ]

  formula_str <- paste(outcome, "~", covariate_formula)
  model <- lm(as.formula(formula_str), data = complete_data)

  residuals_vector <- rep(NA, nrow(data))
  residuals_vector[complete_rows] <- residuals(model)
  return(residuals_vector)
}

cat("Adjusting age_of_onset for covariates using linear regression...\n")
cases$age_of_onset_resid <- get_linear_residuals("age_of_onset", cases)

out_path <- file.path(data_dir, 'phenotypes', 'ad_ageofonset_phenotype_residuals.txt')
write.table(cases, file = out_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nSaved", nrow(cases), "AD cases to:", out_path, "\n")
print(summary(cases$age_of_onset_resid))
