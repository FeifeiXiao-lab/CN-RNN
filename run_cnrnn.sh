#!/usr/bin/env bash
# =============================================================================
# run_cnrnn.sh — CN-RNN end-to-end WES CNV calling pipeline
#
# Usage:
#   bash run_cnrnn.sh [--start N] [--end N]
#
# Options:
#   --start N   Start from step N (1–4). Default: 1
#   --end   N   Stop after step N (1–4). Default: 4
#
# Steps:
#   1  Calculate GC content and mappability (scripts/02-calculate_gc_mapp.sh)
#   2  BAM → Read count matrix              (scripts/03-bam_to_readcount.R)
#   3  QC, normalization, smooth, segment   (scripts/04-qc_norm_smooth_segment.R)
#   4  LSTM prediction                      (scripts/05-cnrnn_predict.py)
#
# Edit the CONFIGURATION block below before running.
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION — edit these paths for your dataset
# =============================================================================

# --- Reference genome build (GRCh38 or GRCh37) ---
GENOME_BUILD="GRCh38"       # GRCh38 or GRCh37

# --- Input files ---
# Target capture BED file (use one of the files in resources/ or your own)
TARGET_BED="resources/hglft_genome_GRCh38.bed"

# Reference FASTA (must match GENOME_BUILD)
REF_FASTA="/path/to/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# Mappability BigWig (100-mer; see README for download instructions)
MAPP_BW="/path/to/reference/hg38.fa.mappability_100bp.bw"

# Sample design matrix (TSV with header: SampleID, BamFileName, Type)
# Type: C = control/reference, T = test/case
SAMPLE_MATRIX="/path/to/sample_matrix.txt"

# --- Trained models ---
DEL_MODEL="models/model_del_20250726.keras"
DUP_MODEL="models/model_dup_20250726.keras"

# --- Output directories ---
OUTPUT_DIR="output"
GC_MAPP_DIR="${OUTPUT_DIR}/gc_mapp"
RC_DIR="${OUTPUT_DIR}/read_counts"
PREPROC_DIR="${OUTPUT_DIR}/preprocessing"
PRED_DIR="${OUTPUT_DIR}/predictions"

# --- Compute resources ---
THREADS=16       # Threads for featureCounts (Step 2)
R_CORES=12       # Parallel cores for CORRseq segmentation (Step 3)

# --- Pipeline parameters ---
ALPHA="0.05,0.01"    # CORRseq significance levels
MIN_MARKERS=5        # Discard CNV calls with fewer markers
MAX_MARKERS=120      # Auto-accept CNV calls with more markers

# =============================================================================
# Argument Parsing
# =============================================================================

START_STEP=1
END_STEP=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START_STEP="$2"; shift 2 ;;
    --end)   END_STEP="$2";   shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# =============================================================================
# Helpers
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
step_enabled() { [[ "$1" -ge "$START_STEP" && "$1" -le "$END_STEP" ]]; }

log "=== CN-RNN Pipeline ==="
log "Steps: ${START_STEP} → ${END_STEP}"
log "Output root: ${OUTPUT_DIR}"
echo

# =============================================================================
# Step 1 — GC Content and Mappability Annotation
# =============================================================================
if step_enabled 1; then
  log "--- Step 1: Calculate GC Content and Mappability ---"
  mkdir -p "${GC_MAPP_DIR}"

  GC_MAPP_FILE="${GC_MAPP_DIR}/Annotated_Targets_${GENOME_BUILD}.txt"

  bash "${SCRIPT_DIR}/scripts/02-calculate_gc_mapp.sh" \
    "${TARGET_BED}" \
    "${REF_FASTA}" \
    "${MAPP_BW}" \
    "${GC_MAPP_FILE}"

  log "Step 1 complete: ${GC_MAPP_FILE}"
  echo
else
  # When skipping Step 1, expect the file to already exist
  GC_MAPP_FILE="${GC_MAPP_DIR}/Annotated_Targets_${GENOME_BUILD}.txt"
  if [[ ! -f "${GC_MAPP_FILE}" ]]; then
    log "ERROR: Skipping Step 1 but GC/Mapp file not found: ${GC_MAPP_FILE}"
    exit 1
  fi
fi

# =============================================================================
# Step 2 — BAM to Read Counts
# =============================================================================
if step_enabled 2; then
  log "--- Step 2: BAM to Read Count Matrix ---"

  Rscript "${SCRIPT_DIR}/scripts/03-bam_to_readcount.R" \
    --sample_matrix "${SAMPLE_MATRIX}" \
    --bed_file      "${TARGET_BED}" \
    --output_dir    "${RC_DIR}" \
    --threads       "${THREADS}" \
    --paired        TRUE

  log "Step 2 complete: ${RC_DIR}/ReadCount_Matrix.txt"
  echo
else
  if [[ ! -f "${RC_DIR}/ReadCount_Matrix.txt" ]]; then
    log "ERROR: Skipping Step 2 but ReadCount_Matrix.txt not found in ${RC_DIR}"
    exit 1
  fi
fi

# =============================================================================
# Step 3 — QC, Normalization, Smoothing, Segmentation
# =============================================================================
if step_enabled 3; then
  log "--- Step 3: QC, Normalization, Smoothing, Segmentation ---"

  Rscript "${SCRIPT_DIR}/scripts/04-qc_norm_smooth_segment.R" \
    --rc_file       "${RC_DIR}/ReadCount_Matrix.txt" \
    --gc_mapp_file  "${GC_MAPP_FILE}" \
    --sample_matrix "${SAMPLE_MATRIX}" \
    --output_dir    "${PREPROC_DIR}" \
    --n_cores       "${R_CORES}" \
    --alpha         "${ALPHA}"

  log "Step 3 complete: ${PREPROC_DIR}"
  echo
else
  if [[ ! -f "${PREPROC_DIR}/smoothed_L2R.h5" ]]; then
    log "ERROR: Skipping Step 3 but smoothed_L2R.h5 not found in ${PREPROC_DIR}"
    exit 1
  fi
fi

# =============================================================================
# Step 4 — LSTM Prediction
# =============================================================================
if step_enabled 4; then
  log "--- Step 4: LSTM CNV Prediction ---"

  python3 "${SCRIPT_DIR}/scripts/05-cnrnn_predict.py" \
    --h5_file     "${PREPROC_DIR}/smoothed_L2R.h5" \
    --calls_dir   "${PREPROC_DIR}/CORRseq_results" \
    --del_model   "$( [[ "${DEL_MODEL}" = /* ]] && echo "${DEL_MODEL}" || echo "${SCRIPT_DIR}/${DEL_MODEL}" )" \
    --dup_model   "$( [[ "${DUP_MODEL}" = /* ]] && echo "${DUP_MODEL}" || echo "${SCRIPT_DIR}/${DUP_MODEL}" )" \
    --output_dir  "${PRED_DIR}" \
    --min_markers "${MIN_MARKERS}" \
    --max_markers "${MAX_MARKERS}" \
    --alpha       "$(echo ${ALPHA} | cut -d',' -f1)"

  log "Step 4 complete: ${PRED_DIR}/Final_CNV_Predictions.csv"
  echo
fi

log "=== Pipeline complete ==="
