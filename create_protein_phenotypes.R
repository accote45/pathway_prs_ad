############################################################
# Residualize UK Biobank protein data and apply RINT
# Modified from residualize_proteins_with_RINT.R - Beatrice W
############################################################

library(data.table)

#-----------------------------------------------------------
# 1. Load protein and covariate data
#-----------------------------------------------------------
raw <- fread("/sc/arion/projects/paul_oreilly/data/ukb/protein/qc/data/protein_measure_baseline.txt")
cov <- fread("/sc/arion/projects/paul_oreilly/data/ukb/protein/qc/data/protein_all_covar.txt")

# APOE genotype (columns: IID, APOE_genotype, e4_copy, e2_copy)
apoe <- fread("/sc/arion/projects/paul_oreilly/lab/shared/UKB_allSample_APOE_genotype.txt")
setnames(apoe, "IID", "sample_id")
cov <- merge(cov, apoe[, .(sample_id, APOE_genotype)], by = "sample_id", all.x = TRUE)

# Load protein QC summary to filter high-missing proteins
missing <- fread("/sc/arion/projects/paul_oreilly/data/ukb/protein/qc/data/protein_qc_summary.txt")
remove <- missing[missing$missing_frac > 0.2, "protein"]

#-----------------------------------------------------------
# 2. Remove proteins with >20% missing values
#-----------------------------------------------------------
d <- raw[, !names(raw) %in% remove$protein, with = FALSE]
setnames(d, "eid", "sample_id")  # rename sample identifier

keep<-fread("/sc/arion/projects/psychgen/projects/prs/extreme_traits/July_allPop_residFirst/ukb18177-allpop-qc.fam")
d<-d[sample_id %in% keep$V1]

#-----------------------------------------------------------
# 4. Prepare covariates
#-----------------------------------------------------------
cov[, Sex := factor(Sex)]
cov[, Batch_info := factor(Batch_info)]
cov[, AC := factor(AC)]
cov[, genotyping_batch := factor(genotyping_batch)]

# APOE: drop e1-containing genotypes, treat blanks as missing, factor
cov[APOE_genotype == "", APOE_genotype := NA]
cov[grepl("e1", APOE_genotype, fixed = TRUE), APOE_genotype := NA]
cov[, APOE_genotype := factor(APOE_genotype)]

# Covariates for residualization
covars <- c(
  "Age", "Age2",
  "Sex",
  "AgeSex", "Age2Sex",
  "Batch_info",
  "AC",
  "genotyping_batch",
  "TBMS",
  "APOE_genotype",
  paste0("PC", 1:15)
)

# Keep only common sample IDs
common_ids <- intersect(d$sample_id, cov$sample_id)

#-----------------------------------------------------------
# 5. Initialize residual table
#-----------------------------------------------------------
# Table with missing values preserved
resid_table1 <- data.table(sample_id = common_ids)

#-----------------------------------------------------------
# 6. Residualize each protein
#-----------------------------------------------------------
# Only residualize these proteins
proteins <- c("NEFL", "GFAP")
protein_idx <- match(proteins, colnames(d))
if (any(is.na(protein_idx))) {
  stop("Protein(s) not found in data: ",
       paste(proteins[is.na(protein_idx)], collapse = ", "))
}

for (i in protein_idx) {

  dat <- merge(d[, c(1, i), with = FALSE], cov, by = "sample_id", all = FALSE)
  dat <- na.omit(dat)  # remove rows with missing values
  
  setnames(dat, colnames(d)[i], "protein_value")
  
  form <- as.formula(paste("protein_value ~", paste(covars, collapse = " + ")))
  fit <- lm(form, data = dat)
  dat[, residual := resid(fit)]
  
  dat_resid <- dat[, .(sample_id, residual)]
  setnames(dat_resid, "residual", colnames(d)[i])
  
  # Merge into main table
  resid_table1 <- merge(resid_table1, dat_resid, by = "sample_id", all.x = TRUE)
}

#-----------------------------------------------------------
# 7. Save residualized protein table
#-----------------------------------------------------------
#fwrite(resid_table1, "residualized_NPX_withNA.txt", sep = "\t", quote = FALSE, na = "NA", col.names = TRUE, row.names = FALSE)

#-----------------------------------------------------------
# 8. Rank-based Inverse Normal Transformation (RINT) function
#-----------------------------------------------------------
rint_dt <- function(dt) {
  dt2 <- copy(dt)
  num_cols <- setdiff(names(dt2), "sample_id")  # skip sample_id
  
  dt2[, (num_cols) := lapply(.SD, function(x) {
    na_idx <- is.na(x)  # preserve NAs
    x_ranked <- rank(x, ties.method = "average", na.last = "keep")
    qnorm((x_ranked - 0.5) / sum(!na_idx))
  }), .SDcols = num_cols]
  
  return(dt2)
}

#-----------------------------------------------------------
# 9. Apply RINT to residualized table
#-----------------------------------------------------------
resid_table1_rint <- rint_dt(resid_table1)

#-----------------------------------------------------------
# 10. Save RINT-transformed table
#-----------------------------------------------------------
fwrite(resid_table1_rint, "residualized_NPX_withNA_RINT.txt", sep = "\t", quote = FALSE, na = "NA", col.names = TRUE, row.names = FALSE)

