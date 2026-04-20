#!/bin/bash
# Script: calculate_features.sh
# Usage: ./calculate_features.sh <INPUT_BED> <REF_FASTA> <MAPPABILITY_BW> <OUTPUT_FILE>

set -e  # Exit immediately if a command exits with a non-zero status

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <INPUT_BED> <REF_FASTA> <MAPPABILITY_BW> <OUTPUT_FILE>"
    exit 1
fi

INPUT_BED=$1
FASTA_FILE=$2
MAPPABILITY_BW=$3
OUTPUT_FILE=$4

# Check for required tools
command -v bedtools >/dev/null 2>&1 || { echo >&2 "Error: bedtools is not installed."; exit 1; }
command -v bigWigAverageOverBed >/dev/null 2>&1 || { echo >&2 "Error: bigWigAverageOverBed is not installed."; exit 1; }

echo "--- Starting Feature Calculation ---"
echo "Target BED: $INPUT_BED"
echo "Reference:  $FASTA_FILE"
echo "Output:     $OUTPUT_FILE"

# Create a temporary workspace
TEMP_PREFIX="tmp_calc_$(date +%s)"

# Step A: Format BED for UCSC Tools
# Adds a unique identifier "region_N" to the 4th column
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "region_"NR}' "$INPUT_BED" > "${TEMP_PREFIX}_formatted.bed"

# Step B: Calculate Mappability (Mean)
echo "[1/3] Calculating Mappability..."
bigWigAverageOverBed "$MAPPABILITY_BW" "${TEMP_PREFIX}_formatted.bed" "${TEMP_PREFIX}_mapp.txt"
# Extract 'mean0' (Column 6)
awk '{print $6}' "${TEMP_PREFIX}_mapp.txt" > "${TEMP_PREFIX}_mapp_col.txt"

# Step C: Calculate GC Content
echo "[2/3] Calculating GC Content..."
# bedtools nuc outputs: Cols 1-3, %AT, %GC... we extract column 5 (%GC)
bedtools nuc -fi "$FASTA_FILE" -bed "$INPUT_BED" | grep -v "^#" | cut -f5 > "${TEMP_PREFIX}_gc_col.txt"

# Step D: Combine and Save
echo "[3/3] Combining Results..."
echo -e "Chr\tStart\tEnd\tName\tGC\tMappability" > "$OUTPUT_FILE"
paste "${TEMP_PREFIX}_formatted.bed" "${TEMP_PREFIX}_gc_col.txt" "${TEMP_PREFIX}_mapp_col.txt" >> "$OUTPUT_FILE"

# Cleanup
rm "${TEMP_PREFIX}_formatted.bed" "${TEMP_PREFIX}_mapp.txt" "${TEMP_PREFIX}_mapp_col.txt" "${TEMP_PREFIX}_gc_col.txt"

echo "Success! Features saved to: $OUTPUT_FILE"
