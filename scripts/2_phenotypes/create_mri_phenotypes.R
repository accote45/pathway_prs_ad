############################################################
# Create UK Biobank MRI endophenotypes for pathway PRS scan
# ----------------------------------------------------------
# Input : UKB_AD_MRI_FSL.csv  (produced by extract_ukb_hipp_wmh.sql)
# Output: residualized + RINT MRI phenotype table for PRSice
#
# Recipe follows the DOWNSTREAM NOTES in extract_ukb_hipp_wmh.sql
# and mirrors the residualize + RINT pattern in create_protein_phenotypes.R
############################################################

library(data.table)

data_dir <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data"

#-----------------------------------------------------------
# 1. Load MRI extract
#-----------------------------------------------------------
mri <- fread(file.path(data_dir, "raw_extracts", "UKB_AD_MRI_FSL.csv"))

# Standardise the sample identifier used downstream (PRSice --ignore-fid keys
# on IID; the SQL extract already emits FID == IID == sample_id).
setnames(mri, "IID", "sample_id")

#-----------------------------------------------------------
# 2. Sample QC exclusions (logged for transparent N-flow)
#-----------------------------------------------------------
# Standard UKB genetic QC, per the SQL downstream notes:
#   Sex == Genetic_Sex
#   is.na(Sex_Chr_Aneuploidy)
#   is.na(Het_Missing_Outlier)
# NOTE: ancestry (EUR) subsetting is NOT done here -- it is applied at PRS
# analysis time via a sample keep-file, to match the sample used in the other
# analyses. Deliberately no White_British filter.
n0 <- nrow(mri)

log_drop <- function(dt, keep, label) {
  kept <- sum(keep, na.rm = TRUE)
  cat(sprintf("QC | %-28s : %6d -> %6d  (dropped %d)\n",
              label, nrow(dt), kept, nrow(dt) - kept))
  dt[keep %in% TRUE]
}

mri <- log_drop(mri, mri$Sex == mri$Genetic_Sex,      "sex concordance")
mri <- log_drop(mri, is.na(mri$Sex_Chr_Aneuploidy),   "no sex-chr aneuploidy")
mri <- log_drop(mri, is.na(mri$Het_Missing_Outlier),  "not het/missing outlier")
cat(sprintf("QC | %-28s : %6d -> %6d\n", "TOTAL", n0, nrow(mri)))

#-----------------------------------------------------------
# 3. Prepare covariates
#-----------------------------------------------------------
mri[, Age  := as.numeric(Age_Imaging)]
mri[, Age2 := Age^2]
mri[, Sex             := factor(Sex)]
mri[, Genotyping_Batch := factor(Genotyping_Batch)]
mri[, Centre_Imaging   := factor(Centre_Imaging)]
mri[, Head_Size_Scaling := as.numeric(Head_Size_Scaling)]

# Baseline covariate set (head-size term added per phenotype below)
base_covars <- c(
  "Age", "Age2",
  "Sex",
  "Genotyping_Batch",
  "Centre_Imaging",
  paste0("PC", 1:15)
)

#-----------------------------------------------------------
# 4. Phenotype specification
#-----------------------------------------------------------
# For each phenotype:
#   transform  : function applied to the raw value before residualization
#   headsize   : include Head_Size_Scaling as a covariate (volume phenotypes)
# Per SQL notes:
#   Hipp_Mean_Vol : residualize on Head_Size_Scaling, then rank-INT
#   WMH_Vol       : log1p(), then rank-INT (no head-size term)
pheno_spec <- list(
  Hipp_Mean_Vol = list(transform = identity,          headsize = TRUE),
  WMH_Vol       = list(transform = function(x) log1p(x), headsize = FALSE)
)

#-----------------------------------------------------------
# 5. Residualize each phenotype (complete cases per phenotype, NAs preserved)
#-----------------------------------------------------------
resid_table <- data.table(sample_id = mri$sample_id)

for (ph in names(pheno_spec)) {
  spec   <- pheno_spec[[ph]]
  covars <- if (spec$headsize) c(base_covars, "Head_Size_Scaling") else base_covars

  dat <- mri[, c("sample_id", ph, covars), with = FALSE]
  setnames(dat, ph, "y")
  dat[, y := spec$transform(as.numeric(y))]
  dat <- na.omit(dat)                       # complete cases for this phenotype

  form <- as.formula(paste("y ~", paste(covars, collapse = " + ")))
  fit  <- lm(form, data = dat)
  dat[, residual := resid(fit)]

  dat_resid <- dat[, .(sample_id, residual)]
  setnames(dat_resid, "residual", ph)
  resid_table <- merge(resid_table, dat_resid, by = "sample_id", all.x = TRUE)

  cat(sprintf("RESID | %-14s : n = %d\n", ph, nrow(dat)))
}

#-----------------------------------------------------------
# 6. Rank-based Inverse Normal Transformation (RINT)
#    (identical to create_protein_phenotypes.R)
#-----------------------------------------------------------
rint_dt <- function(dt) {
  dt2 <- copy(dt)
  num_cols <- setdiff(names(dt2), "sample_id")   # skip sample_id

  dt2[, (num_cols) := lapply(.SD, function(x) {
    na_idx <- is.na(x)                           # preserve NAs
    x_ranked <- rank(x, ties.method = "average", na.last = "keep")
    qnorm((x_ranked - 0.5) / sum(!na_idx))
  }), .SDcols = num_cols]

  dt2
}

resid_table_rint <- rint_dt(resid_table)

#-----------------------------------------------------------
# 7. Save phenotype table (PRSice format: IID + one column per phenotype)
#-----------------------------------------------------------
setnames(resid_table_rint, "sample_id", "IID")

out_path <- file.path(data_dir, "phenotypes", "mri_phenotypes_resid_RINT.txt")
fwrite(resid_table_rint, out_path,
       sep = "\t", quote = FALSE, na = "NA",
       col.names = TRUE, row.names = FALSE)

cat("Wrote:", out_path, "\n")
cat("Phenotype columns:", paste(setdiff(names(resid_table_rint), "IID"), collapse = ", "), "\n")
