#!/usr/bin/env Rscript
# =============================================================================
# Script: 03-bam_to_readcount.R
# Purpose: Count reads per target region from BAM files using featureCounts
#
# Usage:
#   Rscript scripts/03-bam_to_readcount.R \
#       --sample_matrix  <path/to/sample_matrix.txt> \
#       --bed_file       <path/to/targets.bed> \
#       --output_dir     <path/to/output/> \
#       [--threads 16] \
#       [--paired TRUE]
#
# Inputs:
#   sample_matrix : TSV with header; required columns: SampleID, BamFileName, Type
#   bed_file      : BED file (no header), columns: Chr, Start, End
#
# Outputs (written to output_dir):
#   ReadCount_Matrix.txt  — tab-separated read count matrix (Chr, Start, End, samples...)
#   ReadCount_Matrix.rds  — same data as R object
# =============================================================================

suppressPackageStartupMessages({
  library(Rsubread)
  library(dplyr)
})

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
    sample_matrix = get_arg("--sample_matrix"),
    bed_file      = get_arg("--bed_file"),
    output_dir    = get_arg("--output_dir"),
    threads       = as.integer(get_arg("--threads", "16")),
    paired        = as.logical(get_arg("--paired", "TRUE"))
  )
}

params <- parse_args()

# Validate required arguments
required <- c("sample_matrix", "bed_file", "output_dir")
for (arg in required) {
  if (is.null(params[[arg]])) {
    stop(sprintf("Required argument missing: --%s\n", arg))
  }
}
for (arg in c("sample_matrix", "bed_file")) {
  if (!file.exists(params[[arg]])) {
    stop(sprintf("File not found: %s\n", params[[arg]]))
  }
}

cat("=== CN-RNN | Step 3: BAM to Read Counts ===\n")
cat(sprintf("  Sample matrix : %s\n", params$sample_matrix))
cat(sprintf("  Target BED    : %s\n", params$bed_file))
cat(sprintf("  Output dir    : %s\n", params$output_dir))
cat(sprintf("  Threads       : %d\n", params$threads))
cat(sprintf("  Paired-end    : %s\n", params$paired))
cat("\n")

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
if (!dir.exists(params$output_dir)) {
  dir.create(params$output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Load Sample Matrix
# -----------------------------------------------------------------------------
meta <- read.table(params$sample_matrix, header = TRUE, stringsAsFactors = FALSE)

required_cols <- c("SampleID", "BamFileName")
missing_cols  <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(sprintf("Sample matrix is missing required columns: %s\n",
               paste(missing_cols, collapse = ", ")))
}

bam_files  <- meta$BamFileName
sample_ids <- meta$SampleID

missing_bams <- bam_files[!file.exists(bam_files)]
if (length(missing_bams) > 0) {
  stop(sprintf("The following BAM files were not found:\n%s\n",
               paste(missing_bams, collapse = "\n")))
}

cat(sprintf("Samples loaded: %d\n", nrow(meta)))

# -----------------------------------------------------------------------------
# Build SAF Annotation from BED
# -----------------------------------------------------------------------------
bed_df <- read.table(params$bed_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bed_df)[1:3] <- c("Chr", "Start", "End")

saf_df <- data.frame(
  GeneID = paste0("Region_", seq_len(nrow(bed_df))),
  Chr    = bed_df$Chr,
  Start  = bed_df$Start,
  End    = bed_df$End,
  Strand = "*",
  stringsAsFactors = FALSE
)

cat(sprintf("Target regions: %d\n\n", nrow(saf_df)))

# -----------------------------------------------------------------------------
# Run featureCounts
# -----------------------------------------------------------------------------
cat(sprintf("Running featureCounts on %d BAMs...\n", length(bam_files)))

fc <- featureCounts(
  files                 = bam_files,
  annot.ext             = saf_df,
  isGTFAnnotationFile   = FALSE,
  nthreads              = params$threads,
  isPairedEnd           = params$paired,
  requireBothEndsMapped = params$paired,
  allowMultiOverlap     = TRUE
)

# -----------------------------------------------------------------------------
# Format and Save Output
# -----------------------------------------------------------------------------
counts <- fc$counts
colnames(counts) <- sample_ids

final_df <- cbind(bed_df[, 1:3], as.data.frame(counts))

out_txt <- file.path(params$output_dir, "ReadCount_Matrix.txt")
out_rds <- file.path(params$output_dir, "ReadCount_Matrix.rds")

write.table(final_df, out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(final_df, out_rds)

cat(sprintf("\nDone. Outputs:\n  %s\n  %s\n", out_txt, out_rds))
