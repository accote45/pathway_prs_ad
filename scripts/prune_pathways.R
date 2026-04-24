#!/usr/bin/env Rscript
# prune_pathways.R
# Pathway Redundancy Pruning for Pathway PRS Analysis

options(stringsAsFactors = FALSE)

# ==============================================================================
# PARAMETERS — edit before running
# ==============================================================================

# Named list: source prefix (used in SOURCE_PRIORITY) => path to GMT file.
# Add/remove entries to match your available databases.
INPUT_FILES <- list(
  #"GO_BP"                     = "path/to/GO.gmt",
  "SynGO"                     = "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/synGO/qced_data/SynGO.gmt",
  "MitoCarta"                 = "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/mitocarta/qced_data/MitoCarta.gmt",
  #"MSigDB_C2_CP_REACTOME"     = "path/to/MSigDB_C2_CP_REACTOME.gmt",
  #"MSigDB_C2_CP_KEGG"         = "path/to/MSigDB_C2_CP_KEGG.gmt",
  #"MSigDB_C2_CP_WIKIPATHWAYS" = "path/to/MSigDB_C2_CP_WIKIPATHWAYS.gmt",
  #"MSigDB_C2_CP"              = "path/to/MSigDB_C2_CP_other.gmt",
  #"MGI"                       = "path/to/MGI.gmt",
  "DSigDB"                    = "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/dsigdb/qced_data/DSigDB.gmt"
  #"ClinPGx"                   = "path/to/ClinPGx.gmt",
  #"Coexpression_WGCNA"        = "path/to/Coexpression.gmt"
)
OUTPUT_DIR <- "output"

# Similarity thresholds
THRESH_NEAR_DUP    <- 0.8   # min(containment) > this  => near_duplicate
THRESH_NESTED_MAX  <- 0.8   # max(containment) > this  \
THRESH_NESTED_MIN  <- 0.5   # min(containment) < this  |=> nested
THRESH_SIZE_RATIO  <- 0.1   # size_ratio       > this  /

# Source priority: higher number = preferred to keep when resolving near-duplicates.
# Sources are extracted as the prefix of pathway_id before the first "__".
# "Coexpression" is matched by the prefix "Coexpression_".
SOURCE_PRIORITY <- c(
  "SynGO"                     = 10L,
  "MitoCarta"                  = 10L,
  "MSigDB_C2_CP_REACTOME"      = 9L,
  "MSigDB_C2_CP_KEGG"          = 8L,
  "MSigDB_C2_CP_WIKIPATHWAYS"  = 7L,
  "MSigDB_C2_CP"               = 6L,
  "GO_BP"                      = 5L,
  "MGI"                        = 4L,
  "DSigDB"                     = 3L,
  "ClinPGx"                    = 2L,
  "Coexpression"               = 1L   # matched by "Coexpression_" prefix
)

# Parse optional command-line overrides: Rscript prune_pathways.R <input> <outdir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) OUTPUT_DIR <- args[1]

# ==============================================================================
# Helper: look up source priority for a named character vector of sources
# ==============================================================================
get_source_priority <- function(source_vec) {
  vapply(source_vec, function(s) {
    if (s %in% names(SOURCE_PRIORITY)) return(SOURCE_PRIORITY[[s]])
    if (grepl("^Coexpression_", s))    return(SOURCE_PRIORITY[["Coexpression"]])
    0L
  }, integer(1L))
}

# ==============================================================================
# MAIN
# ==============================================================================
main <- function() {

  library(Matrix)
  library(data.table)

  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

  # ============================================================================
  # Step 1: Load and validate
  # ============================================================================
  message("[Step 1] Loading and aggregating input GMT files...")

  read_gmt <- function(path) {
    raw_lines <- readLines(path)
    raw_lines <- raw_lines[nchar(trimws(raw_lines)) > 0L]
    rows <- lapply(raw_lines, function(ln) {
      parts <- strsplit(ln, "\t", fixed = TRUE)[[1L]]
      if (length(parts) < 3L) return(NULL)
      genes <- parts[seq.int(3L, length(parts))]
      genes <- genes[nchar(genes) > 0L]
      if (length(genes) == 0L) return(NULL)
      data.table(pathway_id = parts[1L], placeholder = parts[2L], gene = genes)
    })
    rbindlist(Filter(Negate(is.null), rows))
  }

  all_chunks <- lapply(names(INPUT_FILES), function(src) {
    fpath <- INPUT_FILES[[src]]
    if (!file.exists(fpath)) {
      warning("File not found, skipping: ", fpath)
      return(NULL)
    }
    message("  Reading [", src, "] from: ", fpath)
    chunk <- read_gmt(fpath)
    # Prepend source prefix — avoids collisions across databases
    chunk[, pathway_id := paste0(src, "__", pathway_id)]
    chunk
  })

  dt <- rbindlist(Filter(Negate(is.null), all_chunks))
  message("  Total rows before deduplication: ", nrow(dt))

  # Drop duplicate (pathway_id, gene) pairs
  n_before <- nrow(dt)
  dt <- unique(dt, by = c("pathway_id", "gene"))
  n_dropped <- n_before - nrow(dt)
  if (n_dropped > 0L)
    message("  Dropped ", n_dropped, " duplicate (pathway_id, gene) rows.")

  # Pathway sizes
  psizes <- dt[, .(n_genes = .N), by = pathway_id]

  # Guard: warn and remove pathways outside 10-2000 range
  invalid <- psizes[n_genes < 10L | n_genes > 2000L]
  if (nrow(invalid) > 0L) {
    warning("Removing ", nrow(invalid), " pathway(s) outside size range [10, 2000]: ",
            paste(head(invalid$pathway_id, 5L), collapse = ", "),
            if (nrow(invalid) > 5L) " ..." else "")
    dt     <- dt[!pathway_id %in% invalid$pathway_id]
    psizes <- psizes[!pathway_id %in% invalid$pathway_id]
  }

  # Infer source: everything before the first "__" in pathway_id
  psizes[, source := sub("__.*", "", pathway_id)]

  # One placeholder value per pathway (first occurrence)
  placeholders <- dt[, .(placeholder = placeholder[1L]), by = pathway_id]
  pathway_meta <- merge(psizes, placeholders, by = "pathway_id")

  pathway_ids <- pathway_meta$pathway_id
  N           <- length(pathway_ids)

  message("  ", N, " pathways loaded after validation.")

  # ============================================================================
  # Step 2: Build sparse gene-membership matrix (N × G)
  # ============================================================================
  message("[Step 2] Building sparse gene-membership matrix...")

  all_genes <- unique(dt$gene)
  G         <- length(all_genes)
  gene_idx  <- setNames(seq_len(G),  all_genes)
  pw_idx    <- setNames(seq_len(N),  pathway_ids)

  row_i <- pw_idx[dt$pathway_id]
  col_j <- gene_idx[dt$gene]

  M <- sparseMatrix(
    i        = row_i,
    j        = col_j,
    x        = 1,
    dims     = c(N, G),
    dimnames = list(pathway_ids, all_genes)
  )

  sizes <- setNames(as.integer(rowSums(M)), pathway_ids)
  message("  Matrix: ", N, " pathways × ", G, " genes")

  # ============================================================================
  # Step 3: Compute pairwise intersections via tcrossprod
  # ============================================================================
  message("[Step 3] Computing pairwise intersections...")

  I_mat <- tcrossprod(M)   # N × N sparse symmetric; [a,b] = |genes(a) ∩ genes(b)|

  # ============================================================================
  # Step 4: Compute similarity metrics (upper triangle, nonzero entries only)
  # ============================================================================
  message("[Step 4] Computing similarity metrics...")

  # triu(k=1) excludes the diagonal, so no need to zero it first
  I_upper <- triu(I_mat, k = 1L)
  sp      <- as.data.table(summary(I_upper))   # columns: i, j, x
  setnames(sp, c("i", "j", "intersect"))
  sp[, intersect := as.numeric(intersect)]

  sp[, n_a := sizes[i]]
  sp[, n_b := sizes[j]]
  sp[, containment_a_in_b := intersect / n_a]
  sp[, containment_b_in_a := intersect / n_b]
  sp[, jaccard             := intersect / (n_a + n_b - intersect)]
  sp[, size_ratio          := pmin(n_a, n_b) / pmax(n_a, n_b)]
  sp[, pathway_a           := pathway_ids[i]]
  sp[, pathway_b           := pathway_ids[j]]

  message("  ", nrow(sp), " pairs with nonzero overlap evaluated.")

  # ============================================================================
  # Step 5: Classify pairs (near_duplicate > nested > unflagged)
  # ============================================================================
  message("[Step 5] Classifying pairs...")

  sp[, min_cont := pmin(containment_a_in_b, containment_b_in_a)]
  sp[, max_cont := pmax(containment_a_in_b, containment_b_in_a)]

  sp[, flag_type := fifelse(
    min_cont > THRESH_NEAR_DUP,
    "near_duplicate",
    fifelse(
      max_cont > THRESH_NESTED_MAX &
        min_cont < THRESH_NESTED_MIN &
        size_ratio > THRESH_SIZE_RATIO,
      "nested",
      NA_character_
    )
  )]

  flagged    <- sp[!is.na(flag_type)]
  n_near_dup <- sum(flagged$flag_type == "near_duplicate")
  n_nested   <- sum(flagged$flag_type == "nested")
  message("  Near-duplicate pairs flagged : ", n_near_dup)
  message("  Nested pairs flagged         : ", n_nested)

  # ============================================================================
  # Step 6: Apply greedy single-pass pruning
  # ============================================================================
  message("[Step 6] Applying greedy pruning...")

  status         <- setNames(rep("retained", N), pathway_ids)
  represented_by <- setNames(pathway_ids,         pathway_ids)
  reason         <- setNames(rep("kept",     N), pathway_ids)

  # Precompute source priority for every pathway
  src_prio <- get_source_priority(
    setNames(pathway_meta$source, pathway_meta$pathway_id)
  )

  # --- Near-duplicate pairs: source priority > smaller size > lex id ---
  nd_pairs <- flagged[flag_type == "near_duplicate"]

  for (k in seq_len(nrow(nd_pairs))) {
    pa <- nd_pairs$pathway_a[k]
    pb <- nd_pairs$pathway_b[k]
    if (status[pa] != "retained" || status[pb] != "retained") next

    if      (src_prio[pa] >  src_prio[pb])    { keep <- pa; drop_pw <- pb
    } else if (src_prio[pb] >  src_prio[pa])  { keep <- pb; drop_pw <- pa
    } else if (sizes[pa]   <  sizes[pb])       { keep <- pa; drop_pw <- pb
    } else if (sizes[pb]   <  sizes[pa])       { keep <- pb; drop_pw <- pa
    } else if (pa <= pb)                        { keep <- pa; drop_pw <- pb
    } else                                      { keep <- pb; drop_pw <- pa }

    status[drop_pw]         <- "pruned"
    represented_by[drop_pw] <- keep
    reason[drop_pw]         <- paste0("near_duplicate_of_", keep)
  }

  # --- Nested pairs: always keep the smaller (more specific) pathway ---
  # Sort by descending size of the larger pathway so that broad containment
  # relationships are resolved before narrower ones, preventing chain gaps.
  nest_pairs <- flagged[flag_type == "nested"]
  nest_pairs[, .max_size := pmax(n_a, n_b)]
  setorder(nest_pairs, -.max_size)
  nest_pairs[, .max_size := NULL]

  for (k in seq_len(nrow(nest_pairs))) {
    pa <- nest_pairs$pathway_a[k]
    pb <- nest_pairs$pathway_b[k]
    if (status[pa] != "retained" || status[pb] != "retained") next

    if (sizes[pa] >= sizes[pb]) { keep <- pb; drop_pw <- pa
    } else                       { keep <- pa; drop_pw <- pb }

    status[drop_pw]         <- "pruned"
    represented_by[drop_pw] <- keep
    reason[drop_pw]         <- paste0("nested_in_", keep)
  }

  # Resolve chains: a pathway may have been pruned in favour of another pathway
  # that was itself later pruned.  Follow each represented_by pointer until it
  # lands on a retained pathway.
  repeat {
    pruned_reps <- status[represented_by] == "pruned"
    if (!any(pruned_reps)) break
    represented_by[pruned_reps] <- represented_by[represented_by[pruned_reps]]
  }

  retained_ids <- pathway_ids[status == "retained"]
  n_retained   <- length(retained_ids)
  message("  Retained : ", n_retained, " / ", N,
          "  (", round(100 * (N - n_retained) / N, 1L), "% reduction)")

  # ============================================================================
  # Sanity checks
  # ============================================================================
  message("Running sanity checks...")

  stopifnot("Pruning map row count != N" = length(status) == N)

  bad_rep <- setdiff(unique(represented_by), retained_ids)
  if (length(bad_rep) > 0L)
    stop("represented_by references non-retained pathway id(s): ",
         paste(bad_rep, collapse = ", "))

  if (n_retained >= 2L) {
    M_ret    <- M[retained_ids, , drop = FALSE]
    sz_ret   <- sizes[retained_ids]
    I_ret_up <- triu(tcrossprod(M_ret), k = 1L)
    sr2      <- as.data.table(summary(I_ret_up))
    if (nrow(sr2) > 0L) {
      setnames(sr2, c("i", "j", "isect"))
      sr2[, isect  := as.numeric(isect)]
      sr2[, n_a    := sz_ret[i]]
      sr2[, n_b    := sz_ret[j]]
      sr2[, c_ab   := isect / n_a]
      sr2[, c_ba   := isect / n_b]
      sr2[, min_c  := pmin(c_ab, c_ba)]
      sr2[, max_c  := pmax(c_ab, c_ba)]
      sr2[, sz_rat := pmin(n_a, n_b) / pmax(n_a, n_b)]

      nd_fail   <- sr2[min_c > THRESH_NEAR_DUP]
      nest_fail <- sr2[max_c > THRESH_NESTED_MAX &
                         min_c < THRESH_NESTED_MIN &
                         sz_rat > THRESH_SIZE_RATIO]

      if (nrow(nd_fail) > 0L)
        warning("SANITY CHECK: ", nrow(nd_fail),
                " retained pair(s) still satisfy the near-duplicate criterion.")
      if (nrow(nest_fail) > 0L)
        warning("SANITY CHECK: ", nrow(nest_fail),
                " retained pair(s) still satisfy the nested criterion.")
    }
  }

  message("  Sanity checks complete.")

  # ============================================================================
  # Step 7: Write outputs
  # ============================================================================
  message("[Step 7] Writing outputs to: ", OUTPUT_DIR)

  # 1. pruning_map.tsv
  pruning_map <- data.table(
    pathway_id     = pathway_ids,
    pathway_name   = pathway_ids,   # GMT col 1 serves as both id and name
    source         = pathway_meta$source[match(pathway_ids, pathway_meta$pathway_id)],
    n_genes        = sizes[pathway_ids],
    status         = status[pathway_ids],
    represented_by = represented_by[pathway_ids],
    reason         = reason[pathway_ids]
  )
  fwrite(pruning_map,
         file.path(OUTPUT_DIR, "pruning_map.tsv"),
         sep = "\t", quote = FALSE)

  # 2. pathways_retained.tsv
  dt_ret <- dt[pathway_id %in% retained_ids]
  fwrite(dt_ret,
         file.path(OUTPUT_DIR, "pathways_retained.tsv"),
         sep = "\t", quote = FALSE, col.names = TRUE)

  # 3. pathways_retained.gmt
  ph_lookup    <- setNames(placeholders$placeholder, placeholders$pathway_id)
  ret_by_pw    <- dt_ret[, .(genes = list(gene)), by = pathway_id]
  gmt_lines    <- ret_by_pw[, {
    paste(c(pathway_id, ph_lookup[pathway_id], unlist(genes)), collapse = "\t")
  }, by = pathway_id]$V1
  writeLines(gmt_lines, file.path(OUTPUT_DIR, "pathways_retained.gmt"))

  # 4. redundancy_pairs.tsv
  out_pairs <- flagged[, .(
    pathway_a, pathway_b,
    n_a, n_b, intersect,
    containment_a_in_b, containment_b_in_a,
    jaccard, size_ratio, flag_type
  )]
  fwrite(out_pairs,
         file.path(OUTPUT_DIR, "redundancy_pairs.tsv"),
         sep = "\t", quote = FALSE)

  # 5. summary.txt
  src_in  <- pathway_meta[, .N, by = source][order(-N)]
  src_out <- pathway_meta[pathway_id %in% retained_ids, .N, by = source][order(-N)]

  summary_lines <- c(
    "=== Pathway Pruning Summary ===",
    sprintf("Total input pathways : %d", N),
    "",
    "Input counts by source:",
    sprintf("  %-40s %d", src_in$source, src_in$N),
    "",
    sprintf("Near-duplicate pairs flagged : %d", n_near_dup),
    sprintf("Nested pairs flagged         : %d", n_nested),
    "",
    sprintf("Retained pathways : %d  (%.1f%% reduction from %d)",
            n_retained, 100 * (N - n_retained) / N, N),
    "",
    "Retained counts by source:",
    sprintf("  %-40s %d", src_out$source, src_out$N)
  )
  writeLines(summary_lines, file.path(OUTPUT_DIR, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  message("Done.")
}

if (!interactive()) main()