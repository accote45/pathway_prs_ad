suppressPackageStartupMessages({
  library(MungeSumstats)
  library(data.table)
})

# ── CONFIG ────────────────────────────────────────────────────────────────────
gwas_file      <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/MVP_R4.1000G_AGR.GIA.PheCodes_EndocrineMetabolic_batch1/MVP_R4.1000G_AGR.GIA.PheCodes_EndocrineMetabolic_batch1/MVP_R4.1000G_AGR.Phe_250_2.AFR.GIA.dbGaP.txt.gz"          # path to GWAS summary stats file
output_file    <- "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/t2d_hg19_mvp.txt.gz"          # path for harmonized output file
build_type     <- "lifted"  # "original" or "lifted"
original_build <- "38"        # "37" or "38"

rsid_col   <- "SNP_ID"       # rsID column name; leave "" if absent
chr_col    <- "chrom"        # chromosome column name
pos_col    <- "pos"          # base-pair position column name
A1_col     <- "ea"           # effect allele column name (in derived data)
A2_col     <- "non_ea"       # non-effect allele column name (derived below)
or_col     <- "or"           # odds ratio column name
se_col     <- ""             # standard error column name; leave "" if absent
pval_col   <- "pval"         # p-value column name
eaf_col    <- "af"           # effect allele frequency column name; NULL if absent
n_col      <- "num_samples"  # sample size column name; NULL if absent
info_col   <- NULL           # imputation INFO column name; NULL if absent
# ─────────────────────────────────────────────────────────────────────────────

# ── PRE-PROCESS: derive non-effect allele ─────────────────────────────────────
# The MVP file has ref, alt, and ea (effect allele) columns separately.
# non_ea is whichever of ref/alt is NOT the effect allele.
sumstats_dt <- fread(gwas_file)
sumstats_dt[, non_ea := ifelse(ea == ref, alt, ref)]

# Sanity checks
stopifnot(all(sumstats_dt$ea    %in% c(sumstats_dt$ref, sumstats_dt$alt) |
              sumstats_dt$ea    == sumstats_dt$ref |
              sumstats_dt$ea    == sumstats_dt$alt))
cat("✓ non_ea derived:", nrow(sumstats_dt), "variants\n")
# ─────────────────────────────────────────────────────────────────────────────

# Determine if rsIDs are present
has_rsid <- nzchar(rsid_col)

# Column mapping
# A1 and A2 are swapped to match MungeSumstats expectations
# (MungeSumstats expects A1 = non-effect allele, A2 = effect allele)
col_map <- list(
  SNP  = if (has_rsid) rsid_col else NULL,
  CHR  = chr_col,
  BP   = pos_col,
  A1   = A2_col,   # non-effect allele → A1
  A2   = A1_col,   # effect allele     → A2
  OR   = or_col,
  SE   = if (nzchar(se_col)) se_col else NULL,
  P    = pval_col,
  FRQ  = eaf_col,
  N    = n_col,
  INFO = info_col
)
col_map <- col_map[!sapply(col_map, is.null)]

column_mapping <- data.frame(
  Uncorrected = toupper(unname(unlist(col_map))),
  Corrected   = names(col_map),
  stringsAsFactors = FALSE
)

# Determine genome builds
ref_genome     <- if (original_build == "37") "GRCh37" else "GRCh38"
convert_genome <- if (build_type == "lifted") {
  if (original_build == "37") "GRCh38" else "GRCh37"
} else NULL

# Create log directory
log_dir <- "logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Harmonize
format_sumstats(
  path           = sumstats_dt,   # pass pre-processed data.table
  ref_genome     = ref_genome,
  convert_ref_genome = convert_genome,

  impute_beta=TRUE,
  INFO_filter    = 0,
  FRQ_filter     = 0,
  rmv_chr        = c("X", "Y", "MT"),
  bi_allelic_filter = TRUE,
  frq_is_maf = FALSE,

  allele_flip_check = TRUE,
  allele_flip_drop  = TRUE,
  allele_flip_z     = FALSE,
  allele_flip_frq   = TRUE,

  check_dups     = TRUE,
  sort_coordinates = TRUE,
  N_dropNA       = FALSE,

  save_path      = output_file,
  write_vcf      = FALSE,

  log_folder_ind = FALSE,
  log_folder     = log_dir,
  log_mungesumstats_msgs = TRUE,

  convert_small_p = TRUE,
  convert_large_p = TRUE,
  convert_neg_p   = TRUE,

  snp_ids_are_rs_ids = has_rsid,

  mapping_file   = column_mapping,
  force_new      = TRUE
)

cat("✓ Completed", build_type, "harmonization\n")