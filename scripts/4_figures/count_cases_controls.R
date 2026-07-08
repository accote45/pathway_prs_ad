#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Count AD cases / controls in the FINAL target data of each PRSet analysis.
#
# "Final target data" = the samples PRSice actually used in the regression.
# The authoritative per-sample record of this is the PRSice `.best` file: its
# `In_Regression == "Yes"` rows are exactly the samples that survived --keep,
# the founder filter, genotype QC, and had a non-missing phenotype. We tabulate
# AD status over those samples.
#
# AD status source:
#   - AD case/control runs (AD diagnosis / subset / with-APOE) use phenotype
#     files that already carry an `AD` 0/1 column -> authoritative for that run.
#   - Biomarker / MRI / age-of-onset runs have no AD column in their phenotype;
#     for those we cross-reference the master AD/Dementia file to report the AD
#     composition of the analysed sample (cases vs controls vs unknown status).
#
# Output: per-analysis counts (cases, controls, unknown, total) to stdout and to
# a CSV. Cases = AD==1, controls = AD==0.
#
# Deps: data.table.
#
# Usage (on Minerva, where the results/.best and phenotype files live):
#   Rscript scripts/count_cases_controls.R \
#     [--resultsdir <dir>] [--datadir <dir>] [--outdir figures] \
#     [--adstatus <AD_Dementia file>]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

# ---- args (base-R, no optparse) -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) {
  i <- which(args == flag); if (length(i)) args[i + 1] else default
}
resultsdir <- getarg("--resultsdir",
  "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/results")
datadir <- getarg("--datadir",
  "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data")
outdir <- getarg("--outdir", "figures")
# Master AD/Dementia file: fallback AD-status lookup for the non-AD phenotypes.
adstatus_file <- getarg("--adstatus",
  "/sc/arion/projects/paul_oreilly/lab/shared/pheno/AD_Dementia_28March2024.txt")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- analysis registry ----------------------------------------------------
# For each PRSet run: the .best prefix (matches --out) and the AD-status source.
#   ad_file = NA -> use the master AD/Dementia file (adstatus_file).
#   ad_file set  -> that run's own phenotype file carries the AD column.
reg <- data.table(
  label = c("AD diagnosis (full)", "AD diagnosis (subset 1000/1000)",
            "AD diagnosis (with-APOE)", "AD age of onset",
            "GFAP (plasma)", "NEFL (plasma)",
            "Hippocampal volume (MRI)", "WMH volume (MRI)"),
  out = c("ad_case.control_prset_nothreshold_eur",
          "ad_case.control_prset_nothreshold_eur_subset1000",
          "ad_case.control_prset_nothreshold_eur_withapoe",
          "ad_ageofonset_prset_nothreshold_eur",
          "prset_nothreshold_eur_GFAP", "prset_nothreshold_eur_NEFL",
          "prset_nothreshold_eur_Hipp_Mean_Vol", "prset_nothreshold_eur_WMH_Vol"),
  ad_file = c(
    file.path(datadir, "phenotypes", "ad_phenotype_residuals.txt"),
    file.path(datadir, "phenotypes", "ad_phenotype_residuals_subset1000.txt"),
    file.path(datadir, "phenotypes", "ad_phenotype_residuals_withapoe.txt"),
    NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

# ---- helpers --------------------------------------------------------------
# IIDs used in the regression, from the .best file (In_Regression == Yes).
# The PRSet .best has one PRS column per pathway (17k), so for the full sample
# it is tens of GB and cannot be memory-mapped. We therefore (a) read ONLY the
# header line to locate the IID / In_Regression columns, then (b) stream those
# columns with awk, filtering to In_Regression == "Yes" so only a tiny list of
# IIDs is ever loaded into R.
used_iids <- function(best_path) {
  if (!file.exists(best_path)) return(NULL)
  hdr <- tryCatch(strsplit(trimws(readLines(best_path, n = 1)), "\\s+")[[1]],
                  error = function(e) NULL)
  if (is.null(hdr)) return(NULL)
  ic <- match("In_Regression", hdr); ii <- match("IID", hdr)
  if (is.na(ii)) return(NULL)
  cmd <- if (is.na(ic))                     # no In_Regression col: take all IIDs
    sprintf("awk 'NR>1 {print $%d}' %s", ii, shQuote(best_path))
  else
    sprintf("awk 'NR>1 && $%d==\"Yes\" {print $%d}' %s", ic, ii, shQuote(best_path))
  v <- tryCatch(fread(cmd = cmd, header = FALSE)[[1]], error = function(e) NULL)
  if (is.null(v)) return(NULL)
  unique(as.character(v))
}

# sample_id -> AD (0/1) lookup from a phenotype/status file.
ad_lookup <- function(path, id_col = "sample_id", ad_col = "AD") {
  if (!file.exists(path)) stop("AD-status file not found: ", path)
  d <- fread(path)
  if (!all(c(id_col, ad_col) %in% names(d)))
    stop("Expected columns '", id_col, "' and '", ad_col, "' in ", path,
         " (found: ", paste(head(names(d), 20), collapse = ", "), ")")
  out <- data.table(id = as.character(d[[id_col]]), AD = suppressWarnings(as.integer(d[[ad_col]])))
  out[!is.na(id)]
}

master_ad <- ad_lookup(adstatus_file)

# ---- tabulate per analysis ------------------------------------------------
rows <- list()
for (i in seq_len(nrow(reg))) {
  lab  <- reg$label[i]
  best <- file.path(resultsdir, paste0(reg$out[i], ".best"))
  message("[", i, "/", nrow(reg), "] ", lab,
          " -- streaming ", basename(best), " (large files take a few minutes) ...")
  iids <- used_iids(best)
  if (is.null(iids)) {
    message("SKIP (.best missing or unreadable): ", best)
    rows[[lab]] <- data.table(analysis = lab, cases = NA, controls = NA,
                              unknown_status = NA, total_used = NA,
                              ad_source = NA_character_, best_found = FALSE)
    next
  }
  # choose AD-status source
  if (!is.na(reg$ad_file[i]) && file.exists(reg$ad_file[i])) {
    look <- ad_lookup(reg$ad_file[i]); src <- basename(reg$ad_file[i])
  } else {
    look <- master_ad; src <- paste0(basename(adstatus_file), " (master)")
  }
  ad <- look$AD[match(iids, look$id)]
  rows[[lab]] <- data.table(
    analysis = lab,
    cases          = sum(ad == 1, na.rm = TRUE),
    controls       = sum(ad == 0, na.rm = TRUE),
    unknown_status = sum(is.na(ad)),
    total_used     = length(iids),
    ad_source      = src,
    best_found     = TRUE)
}
res <- rbindlist(rows)
res[, pct_case := ifelse(is.na(cases), NA, round(100 * cases / (cases + controls), 1))]

# ---- report ---------------------------------------------------------------
cat("\nAD case/control composition of the final target data per PRSet analysis\n")
cat("(cases = AD==1, controls = AD==0; 'unknown' = In_Regression sample with no",
    "AD-status match)\n\n")
print(res[, .(analysis, cases, controls, unknown_status, total_used, pct_case, ad_source)],
      row.names = FALSE)

out_csv <- file.path(outdir, "prset_case_control_counts.csv")
fwrite(res, out_csv)
cat("\nWrote", out_csv, "\n")

miss <- res[best_found == FALSE, analysis]
if (length(miss))
  cat("\nNote: no .best file for:", paste(miss, collapse = "; "),
      "\n  (re-run those PRSet jobs, or point --resultsdir at the right folder;",
      "\n   PRSice writes <out>.best unless regression was skipped.)\n")
cat("\nReminder: age-of-onset is cases-only by construction, so ~0 controls is",
    "expected there.\n")
