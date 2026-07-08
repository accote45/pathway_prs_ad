# ---------------------------------------------------------------------------
# Subset the AD case/control phenotype to 1000 cases + 1000 controls
#
# Repeats the PRSet AD dx association in a reduced sample. We DO NOT re-run the
# covariate residualization: the AD_resid column from create_ad_phenotype.R was
# fit in the full sample, where the nuisance covariate coefficients (Sex, Age,
# Age2, Batch, APOE_genotype, PC1-15) are estimated far more precisely and
# without the sparse-cell / quasi-separation problems that a 2000-person re-fit
# would introduce. The subset is a random draw, so the full-sample adjustment
# is unbiased for it. We simply subset the rows and keep the existing residual.
#
# Outputs:
#   ad_phenotype_residuals_subset1000.txt  -> --pheno for PRSet
#   ad_subset1000_keep.txt (FID IID)       -> --keep for PRSet
# ---------------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(dplyr)

data_dir <- '/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data'

set.seed(42)               # reproducible subsetting
n_cases    <- 1000
n_controls <- 1000

# ---------------------------------------------------------------------------
# Read the full phenotype file produced by create_ad_phenotype.R. It already
# carries AD status and the full-sample residual AD_resid.
# ---------------------------------------------------------------------------
full_pheno <- fread(file.path(data_dir, 'phenotypes', 'ad_phenotype_residuals.txt'))

# Only individuals with a usable (non-NA) full-sample residual are eligible.
dat <- full_pheno[!is.na(full_pheno$AD_resid), ]

cat("Eligible pool: ",
    sum(dat$AD == 1), " cases / ", sum(dat$AD == 0), " controls\n", sep = "")

if (sum(dat$AD == 1) < n_cases || sum(dat$AD == 0) < n_controls) {
  stop("Not enough cases/controls to draw the requested subset.")
}

# ---------------------------------------------------------------------------
# Draw 1000 cases + 1000 controls (keeping the full-sample AD_resid as-is)
# ---------------------------------------------------------------------------
sampled_cases    <- sample(dat$sample_id[dat$AD == 1], n_cases)
sampled_controls <- sample(dat$sample_id[dat$AD == 0], n_controls)
sampled_ids      <- c(sampled_cases, sampled_controls)

sub <- dat[dat$sample_id %in% sampled_ids, ]
cat("Subset drawn: ",
    sum(sub$AD == 1), " cases / ", sum(sub$AD == 0), " controls\n", sep = "")

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
# Phenotype file (same layout as ad_phenotype_residuals.txt; PRSet reads AD_resid)
out_pheno <- file.path(data_dir, 'phenotypes', 'ad_phenotype_residuals_subset1000.txt')
write.table(sub, file = out_pheno, sep = "\t", row.names = FALSE, quote = FALSE)

# Keep file for --keep (FID IID). UKB fam files use FID == IID.
keep <- data.table(FID = sub$sample_id, IID = sub$sample_id)
out_keep <- file.path(data_dir, 'samples', 'ad_subset1000_keep.txt')
write.table(keep, file = out_keep, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

cat("\nWrote:\n  ", out_pheno, "\n  ", out_keep, "\n", sep = "")
cat("Subset AD residual summary:\n")
print(summary(sub$AD_resid))
