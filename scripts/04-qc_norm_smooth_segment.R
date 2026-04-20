#!/usr/bin/env Rscript
# =============================================================================
# Script: 04-qc_norm_smooth_segment.R
# Purpose: QC, GC/mappability normalization, smoothing, and CORRseq segmentation
#
# Usage:
#   Rscript scripts/04-qc_norm_smooth_segment.R \
#       --rc_file       <ReadCount_Matrix.txt> \
#       --gc_mapp_file  <Annotated_Targets.txt> \
#       --sample_matrix <sample_matrix.txt> \
#       --output_dir    <path/to/output/> \
#       [--n_cores 12] \
#       [--alpha "0.05,0.01"]
#
# Inputs:
#   rc_file       : Read count matrix from Step 3 (Chr, Start, End, samples...)
#   gc_mapp_file  : GC + mappability annotation file from Step 2
#   sample_matrix : TSV with columns SampleID, BamFileName, Type (C=control, T=case)
#
# Outputs (written to output_dir):
#   smoothed_L2R.h5          — HDF5 with smoothed log2 ratios + reference QC
#   smoothed_L2R.rds         — same as RDS backup
#   CORRseq_results/         — per-chromosome segmentation CSVs
#     chr[1-22].alpha.[0.05|0.01].csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(rhdf5)
  library(GenomicRanges)
  library(Rcpp)
  library(parallel)
})

# -----------------------------------------------------------------------------
# Locate Script Directory for Relative Library Paths
# -----------------------------------------------------------------------------
args_full   <- commandArgs(trailingOnly = FALSE)
file_flag   <- args_full[grep("--file=", args_full)]
script_dir  <- if (length(file_flag) > 0) {
  dirname(normalizePath(sub("--file=", "", file_flag)))
} else {
  file.path(getwd(), "scripts")   # fallback when sourced interactively
}

lib_dir     <- file.path(script_dir, "library_functions")
corrseq_dir <- file.path(lib_dir, "CORRseq_functions")
modsara_dir <- file.path(lib_dir, "modSaRa_vendor")

# -----------------------------------------------------------------------------
# Load Library Functions
# -----------------------------------------------------------------------------
cat("Loading library functions...\n")

for (f in c("library_functions.R", "Library_size_adjust.R")) {
  path <- file.path(lib_dir, f)
  if (!file.exists(path)) stop(sprintf("Library file not found: %s", path))
  source(path)
}

# Load vendored modSaRa R functions from the repo. We intentionally skip the
# package-generated wrappers that expect an installed package namespace and use
# the local smooth.cpp implementation below instead.
modsara_files <- c(
  "CNVcluster.R",
  "CNVout.R",
  "QuanNorm.R",
  "SARA.R",
  "SARAp.R",
  "estimateSigma.R",
  "fInverse.R",
  "gausianMixure.R",
  "getOneBIC.R",
  "localDiagnostic.R",
  "localMax.R",
  "lrr.R",
  "map.R",
  "modifiedSaRa.R",
  "multiSaRa.R"
)

for (f in modsara_files) {
  path <- file.path(modsara_dir, f)
  if (!file.exists(path)) stop(sprintf("Vendored modSaRa file not found: %s", path))
  source(path)
}

for (f in c("CNVcluster_new.R", "CNVout.R", "gausianMixture.R",
            "modifiedSaRa.R", "multiSaRa.R", "getOneBIC.R",
            "Normalization.R", "localMax.R", "process_chromosome_CORRseq.R")) {
  path <- file.path(corrseq_dir, f)
  if (!file.exists(path)) stop(sprintf("CORRseq function file not found: %s", path))
  source(path)
}

smooth_cpp <- file.path(lib_dir, "smooth.cpp")
if (!file.exists(smooth_cpp)) stop(sprintf("smooth.cpp not found: %s", smooth_cpp))
sourceCpp(smooth_cpp)

# -----------------------------------------------------------------------------
# Argument Parsing
# -----------------------------------------------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) == 0) return(default)
    if (idx >= length(args)) stop(paste("Flag", flag, "requires a value."))
    args[idx + 1]
  }
  list(
    rc_file       = get_arg("--rc_file"),
    gc_mapp_file  = get_arg("--gc_mapp_file"),
    sample_matrix = get_arg("--sample_matrix"),
    output_dir    = get_arg("--output_dir"),
    n_cores       = as.integer(get_arg("--n_cores", "12")),
    alpha         = get_arg("--alpha", "0.05,0.01")
  )
}

params <- parse_args()

required <- c("rc_file", "gc_mapp_file", "sample_matrix", "output_dir")
for (arg in required) {
  if (is.null(params[[arg]])) stop(sprintf("Required argument missing: --%s\n", arg))
}
for (arg in c("rc_file", "gc_mapp_file", "sample_matrix")) {
  if (!file.exists(params[[arg]])) stop(sprintf("File not found: %s\n", params[[arg]]))
}

alpha_vector <- as.numeric(strsplit(params$alpha, ",")[[1]])
n_cores      <- min(params$n_cores, max(1L, detectCores() - 1L))

cat("=== CN-RNN | Step 4: QC, Normalization, Smoothing, Segmentation ===\n")
cat(sprintf("  Read count file  : %s\n", params$rc_file))
cat(sprintf("  GC/Mapp file     : %s\n", params$gc_mapp_file))
cat(sprintf("  Sample matrix    : %s\n", params$sample_matrix))
cat(sprintf("  Output dir       : %s\n", params$output_dir))
cat(sprintf("  Cores            : %d\n", n_cores))
cat(sprintf("  Alpha levels     : %s\n", params$alpha))
cat("\n")

if (!dir.exists(params$output_dir)) dir.create(params$output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
cat("Loading read count matrix...\n")
counts_dt <- fread(params$rc_file)

cat("Loading GC/Mappability annotation...\n")
targets_dt <- fread(params$gc_mapp_file)

# Alignment check
if (nrow(counts_dt) != nrow(targets_dt)) {
  stop(sprintf(
    "Row count mismatch: read count matrix (%d rows) vs targets file (%d rows).\n",
    nrow(counts_dt), nrow(targets_dt)
  ))
}
if (!all(counts_dt$Chr == targets_dt$Chr)) {
  stop("Chromosome mismatch between read count matrix and targets file. Ensure files are sorted identically.\n")
}
if (!all(counts_dt$Start == targets_dt$Start)) {
  stop("Start coordinate mismatch between read count matrix and targets file. Ensure files are sorted identically.\n")
}
if (!all(counts_dt$End == targets_dt$End)) {
  stop("End coordinate mismatch between read count matrix and targets file. Ensure files are sorted identically.\n")
}

# -----------------------------------------------------------------------------
# Build Count Matrix and GRanges Reference
# -----------------------------------------------------------------------------
sample_cols <- 4:ncol(counts_dt)
read_counts <- as.matrix(counts_dt[, ..sample_cols])
rownames(read_counts) <- paste0("region_", seq_len(nrow(read_counts)))

# Support both capitalisation styles for column names
gc_col   <- intersect(c("GC",   "gc"),           colnames(targets_dt))[1]
mapp_col <- intersect(c("Mappability", "mapp"),   colnames(targets_dt))[1]
if (is.na(gc_col))   stop("GC column not found in gc_mapp_file (expected 'GC' or 'gc').\n")
if (is.na(mapp_col)) stop("Mappability column not found in gc_mapp_file (expected 'Mappability' or 'mapp').\n")

GC_Mapp <- GRanges(
  seqnames = targets_dt$Chr,
  ranges   = IRanges(start = targets_dt$Start, end = targets_dt$End),
  gc       = targets_dt[[gc_col]],
  mapp     = targets_dt[[mapp_col]]
)

# -----------------------------------------------------------------------------
# QC Filter
# -----------------------------------------------------------------------------
cat("Running QC filters...\n")
after_qc <- qc_exons(
  Y                     = read_counts,
  ref                   = GC_Mapp,
  cov_threshold         = c(20, 2000),
  length_threshold      = c(20, 19999),
  mappability_threshold = 0.9,
  gc_threshold          = c(0.2, 0.8)
)
cat(sprintf("Regions retained after QC: %d / %d\n\n", nrow(after_qc$Y_qc), nrow(read_counts)))

# -----------------------------------------------------------------------------
# Library Size Adjustment and EMRC
# -----------------------------------------------------------------------------
cat("Library size adjustment...\n")
adjusted_Y_qc <- Library_size_adjust_by_median_Read_depth(read_counts = after_qc$Y_qc)
EMRC <- adjusted_Y_qc / width(after_qc$ref_qc)

# -----------------------------------------------------------------------------
# Identify Reference vs. Case Samples
# -----------------------------------------------------------------------------
cat("Identifying reference and case samples...\n")
sample_df            <- fread(params$sample_matrix, header = TRUE)
available_samples    <- colnames(EMRC)
reference_sample_ids <- intersect(sample_df$SampleID[sample_df$Type == "C"], available_samples)
case_sample_ids      <- intersect(sample_df$SampleID[sample_df$Type == "T"], available_samples)

cat(sprintf("  Reference (C): %d\n", length(reference_sample_ids)))
cat(sprintf("  Case      (T): %d\n\n", length(case_sample_ids)))

if (length(reference_sample_ids) < 5) {
  warning("Fewer than 5 reference samples found — normalization may be unstable.")
}

# -----------------------------------------------------------------------------
# GC/Mappability/Size Normalization → Log2 Ratio
# -----------------------------------------------------------------------------
cat("Running normalization...\n")
L2R <- normalization_bin_ref(
  Y        = EMRC,
  bin_size = width(after_qc$ref_qc),
  gc       = after_qc$ref_qc$gc,
  map      = after_qc$ref_qc$mapp,
  ref_id   = reference_sample_ids
)

# -----------------------------------------------------------------------------
# Smoothing (C++ Rcpp, window R=5)
# -----------------------------------------------------------------------------
cat("Smoothing log2 ratios (window R=5)...\n")
smoothed_L2R <- smooth(lrr = L2R$log2Rdata, R = 5)

# -----------------------------------------------------------------------------
# Save HDF5 and RDS
# -----------------------------------------------------------------------------
h5_path  <- file.path(params$output_dir, "smoothed_L2R.h5")
rds_path <- file.path(params$output_dir, "smoothed_L2R.rds")

cat(sprintf("Saving HDF5: %s\n", h5_path))
ref_qc_df          <- as.data.frame(after_qc$ref_qc)
ref_qc_df$seqnames <- as.character(ref_qc_df$seqnames)
ref_qc_df$strand   <- as.character(ref_qc_df$strand)

if (file.exists(h5_path)) file.remove(h5_path)
h5createFile(h5_path)
h5write(smoothed_L2R,         h5_path, "smoothed_L2R")
h5write(ref_qc_df,            h5_path, "ref_qc")
h5write(case_sample_ids,      h5_path, "case_sample_ids")
h5write(reference_sample_ids, h5_path, "reference_sample_ids")
H5close()

saveRDS(list(
  smoothed_L2R = smoothed_L2R,
  ref_qc       = after_qc$ref_qc,
  case_ids     = case_sample_ids,
  ref_ids      = reference_sample_ids
), rds_path)

cat("Pre-processing complete.\n\n")

# -----------------------------------------------------------------------------
# CORRseq Segmentation (Parallel by Chromosome)
# -----------------------------------------------------------------------------
seg_dir <- file.path(params$output_dir, "CORRseq_results")
if (!dir.exists(seg_dir)) dir.create(seg_dir, recursive = TRUE)

ref_qc_df <- as.data.frame(after_qc$ref_qc)

first_chr   <- as.character(ref_qc_df$seqnames[1])
chromosomes <- if (grepl("^chr", first_chr)) paste0("chr", 1:22) else as.character(1:22)

cat(sprintf("Running CORRseq segmentation on %d cores...\n", n_cores))
cat(sprintf("Output: %s\n\n", seg_dir))

results <- mclapply(
  chromosomes,
  function(chr) {
    tryCatch({
      process_chromosome_CORRseq(
        chr          = chr,
        log2R        = smoothed_L2R,
        ref_qc       = ref_qc_df,
        output_dir   = seg_dir,
        alpha_vector = alpha_vector
      )
      sprintf("%s : OK", chr)
    }, error = function(e) {
      sprintf("%s : FAILED — %s", chr, e$message)
    })
  },
  mc.cores = n_cores
)

cat("--- Segmentation Summary ---\n")
invisible(lapply(results, function(r) cat(sprintf("  %s\n", r))))
cat(sprintf("\nStep 4 complete. Results in: %s\n", seg_dir))
