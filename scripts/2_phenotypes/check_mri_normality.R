############################################################
# Diagnostic: does each MRI phenotype actually need RINT?
# ----------------------------------------------------------
# Residualizes exactly as create_mri_phenotypes.R does, then
# reports distribution shape of the RESIDUALS (pre-RINT).
# Decide per-phenotype -- do NOT assume proteins' RINT carries over.
#
# Rule of thumb (residuals):
#   |skew| < 0.5 and |excess kurtosis| < 1  -> roughly normal, RINT optional
#   |skew| 0.5-1                            -> mild, RINT defensible
#   |skew| > 1                              -> skewed, transform/RINT warranted
# NOTE: ignore Shapiro-Wilk p-values at this N -- they reject trivial
#       deviations. Trust skew/kurtosis + the Q-Q plot.
############################################################

library(data.table)

data_dir <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data"

mri <- fread(file.path(data_dir, "UKB_AD_MRI_FSL.csv"))
setnames(mri, "IID", "sample_id")

# --- same QC as the phenotype script ---
mri <- mri[mri$Sex == mri$Genetic_Sex]
mri <- mri[is.na(mri$Sex_Chr_Aneuploidy)]
mri <- mri[is.na(mri$Het_Missing_Outlier)]

# --- covariates ---
mri[, Age  := as.numeric(Age_Imaging)]
mri[, Age2 := Age^2]
mri[, Sex              := factor(Sex)]
mri[, Genotyping_Batch := factor(Genotyping_Batch)]
mri[, Centre_Imaging   := factor(Centre_Imaging)]
mri[, Head_Size_Scaling := as.numeric(Head_Size_Scaling)]

base_covars <- c("Age", "Age2", "Sex", "Genotyping_Batch",
                 "Centre_Imaging", paste0("PC", 1:15))

# transform applied before residualization (mirror the real script)
pheno_spec <- list(
  Hipp_Mean_Vol_raw = list(col = "Hipp_Mean_Vol", transform = identity,            headsize = TRUE),
  WMH_Vol_raw       = list(col = "WMH_Vol",       transform = identity,            headsize = FALSE),
  WMH_Vol_log1p     = list(col = "WMH_Vol",       transform = function(x) log1p(x), headsize = FALSE)
)

# lightweight skew / excess-kurtosis (no extra packages)
skew <- function(x) { x <- x[!is.na(x)]; m <- mean(x); mean((x-m)^3)/mean((x-m)^2)^1.5 }
ekurt <- function(x) { x <- x[!is.na(x)]; m <- mean(x); mean((x-m)^4)/mean((x-m)^2)^2 - 3 }

# --- audit: how many non-NA values does each covariate/phenotype have? ---
audit_cols <- c(base_covars, "Head_Size_Scaling", "Hipp_Mean_Vol", "WMH_Vol")
cat("\nNon-NA counts (of", nrow(mri), "rows after QC):\n")
for (cn in audit_cols) {
  if (!cn %in% names(mri)) { cat(sprintf("  %-20s MISSING FROM DATA\n", cn)); next }
  nn <- sum(!is.na(mri[[cn]]))
  flag <- if (nn == 0) "  <-- ALL NA" else if (nn < 0.5 * nrow(mri)) "  <-- mostly NA" else ""
  cat(sprintf("  %-20s %8d%s\n", cn, nn, flag))
}
cat("\n")

pdf(file.path(data_dir, "mri_residual_qqplots.pdf"), width = 7, height = 5)

cat(sprintf("%-18s %8s %8s %10s   %s\n",
            "phenotype", "skew", "e.kurt", "n", "verdict"))
cat(strrep("-", 70), "\n")

for (nm in names(pheno_spec)) {
  spec   <- pheno_spec[[nm]]
  covars <- if (spec$headsize) c(base_covars, "Head_Size_Scaling") else base_covars

  dat <- mri[, c("sample_id", spec$col, covars), with = FALSE]
  setnames(dat, spec$col, "y")
  dat[, y := spec$transform(as.numeric(y))]

  # Drop covariates that are entirely NA *before* na.omit, otherwise a single
  # all-NA column (e.g. mis-indexed PC1) deletes every row.
  covars_present <- covars[vapply(covars,
    function(cv) sum(!is.na(dat[[cv]])) > 0, logical(1))]

  dat <- na.omit(dat[, c("sample_id", "y", covars_present), with = FALSE])
  dat <- droplevels(dat)   # shed factor levels absent in this subset

  # drop covariates that don't vary here (single-level factor / constant)
  usable <- covars_present[vapply(covars_present, function(cv) {
    v <- dat[[cv]]
    if (is.factor(v)) nlevels(v) > 1 else length(unique(v)) > 1
  }, logical(1))]
  dropped <- setdiff(covars, usable)
  if (length(dropped))
    cat(sprintf("  [%s] dropped unusable covariate(s): %s\n",
                nm, paste(dropped, collapse = ", ")))

  form <- as.formula(paste("y ~", paste(usable, collapse = " + ")))
  r <- resid(lm(form, data = dat))

  s <- skew(r); k <- ekurt(r)
  verdict <- if (abs(s) < 0.5 && abs(k) < 1) "~normal (RINT optional)"
             else if (abs(s) < 1)            "mild skew (RINT defensible)"
             else                            "skewed (transform/RINT warranted)"

  cat(sprintf("%-18s %8.3f %8.3f %10d   %s\n", nm, s, k, length(r), verdict))

  qqnorm(r, main = sprintf("Q-Q: %s  (skew=%.2f, e.kurt=%.2f)", nm, s, k),
         pch = ".", cex = 2)
  qqline(r, col = "red")
}

dev.off()
cat("\nQ-Q plots written to:", file.path(data_dir, "mri_residual_qqplots.pdf"), "\n")
