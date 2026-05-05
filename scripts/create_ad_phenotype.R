# create master phenotype file
# read in case/control status + age of onset
# read in covariates (age, sex, PCs)
# get phenotype residuals

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


}
















# Adjust outcomes for covariates, obtain residuals

# Define covariates
covariates <- c("Sex", "Age", "Age2", "Batch", "Centre", 
                paste0("PC", 1:15))

covariate_formula <- paste(covariates, collapse = " + ")











# Create formula for covariates
covariate_formula <- paste(covariates, collapse = " + ")

# Function to get residuals from regression
get_residuals <- function(outcome, data) {
  # Create complete cases for this outcome
  complete_data <- data[!is.na(data[[outcome]]), ]
  
  # Create formula
  formula_str <- paste(outcome, "~", covariate_formula)
  
  # Fit model
  model <- lm(as.formula(formula_str), data = complete_data)
  
  # Get residuals and put them back in original data structure
  residuals_vector <- rep(NA, nrow(data))
  residuals_vector[!is.na(data[[outcome]])] <- residuals(model)
  
  return(residuals_vector)
}

# For binary outcomes, use logistic regression
get_logistic_residuals <- function(outcome, data) {
  # Create complete cases for this outcome
  complete_data <- data[!is.na(data[[outcome]]), ]
  
  # Create formula
  formula_str <- paste(outcome, "~", covariate_formula)
  
  # Fit logistic model
  model <- glm(as.formula(formula_str), data = complete_data, family = binomial)
  
  # Get Pearson residuals
  residuals_vector <- rep(NA, nrow(data))
  residuals_vector[!is.na(data[[outcome]])] <- residuals(model, type = "pearson")
  
  return(residuals_vector)
}

# Apply logistic regression adjustment to binary outcomes
    final_pheno[[paste0(outcome, "_resid")]] <- get_logistic_residuals(outcome, final_pheno)

# Save the adjusted phenotype file
write.table(final_pheno, 
            file = '/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb_phenofile_forprset.txt', 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary of adjustments
cat("\nAdjusted outcomes available:\n")
resid_cols <- colnames(final_pheno)[grepl("_resid$", colnames(final_pheno))]
cat(paste(resid_cols, collapse = "\n"))




