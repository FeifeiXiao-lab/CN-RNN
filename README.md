# CN-RNN: Deep Learning-Based CNV Calling for Whole Exome Sequencing

CN-RNN is a WES copy number variant (CNV) calling pipeline that combines classical signal processing with a dual-input LSTM neural network. It takes aligned BAM files as input and produces a table of confirmed CNV calls (deletions and duplications) per sample.

## Pipeline Overview

```
BAM files + Reference genome
   |
   v
[Step 1] GC content + mappability annotation
   |
   v
[Step 2] Read count matrix (featureCounts / Rsubread)
   |
   v
[Step 3] QC -> normalization -> smoothing -> CORRseq segmentation
   |
   v
[Step 4] LSTM refinement -> Final_CNV_Predictions.csv
```

## Repository Structure

```
CN-RNN/
├── run_cnrnn.sh                     # Main pipeline runner
├── cnrnn_full_environment.yml       # Conda environment specification
├── .gitignore
├── resources/
│   ├── hglft_genome_GRCh38.bed     # Exome target BED (GRCh38)
│   └── targets_GRCh37.bed          # Exome target BED (GRCh37)
├── models/
│   ├── model_del_20250726.keras    # Trained deletion LSTM model
│   └── model_dup_20250726.keras    # Trained duplication LSTM model
└── scripts/
    ├── 02-calculate_gc_mapp.sh     # Step 1: GC + mappability annotation
    ├── 03-bam_to_readcount.R       # Step 2: BAM -> read count matrix
    ├── 04-qc_norm_smooth_segment.R # Step 3: QC, normalization, smoothing, segmentation
    ├── 05-cnrnn_predict.py         # Step 4: LSTM prediction
    └── library_functions/          # R library functions (CORRseq + vendored modSaRa)
```

## Requirements

- Linux (tested on CentOS/RHEL 7+, Ubuntu 20.04+)
- Conda or Mamba

### Install the environment

```bash
conda env create -f cnrnn_full_environment.yml
conda activate cnrnn_env
```

This installs all R packages (including Bioconductor), Python packages (TensorFlow 2.15, Keras 2.15), and CLI tools (samtools, bedtools, bigWigAverageOverBed).

The `modSaRa` segmentation library is vendored directly in the repository — no separate installation step is required.

## Quick Start

### 1. Download reference data

You need two reference files matching your genome build:

**GRCh38 (hg38):**

```bash
# Reference FASTA
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

# Mappability BigWig (100-mer)
wget https://zenodo.org/records/5521424/files/hg38.fa.mappability_100bp.bw
```

**GRCh37 (hg19):**

```bash
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz && samtools faidx hs37d5.fa

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
```

### 2. Prepare your sample design matrix

Create a tab-separated file with columns `SampleID`, `BamFileName`, and `Type`:

```
SampleID	BamFileName	Type
CTRL_01	/path/to/CTRL_01.bam	C
CTRL_02	/path/to/CTRL_02.bam	C
CTRL_03	/path/to/CTRL_03.bam	C
CASE_01	/path/to/CASE_01.bam	T
CASE_02	/path/to/CASE_02.bam	T
```

- `Type = C` : control/reference samples (at least 5 recommended)
- `Type = T` : test/case samples (CNVs will be called in these)
- BAM files must be coordinate-sorted and indexed (`.bai` present)

### 3. Configure the pipeline

Edit the **CONFIGURATION** section at the top of `run_cnrnn.sh`:

```bash
GENOME_BUILD="GRCh38"
TARGET_BED="resources/hglft_genome_GRCh38.bed"
REF_FASTA="/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa"
MAPP_BW="/path/to/hg38.fa.mappability_100bp.bw"
SAMPLE_MATRIX="/path/to/sample_matrix.txt"
```

### 4. Run

```bash
conda activate cnrnn_env
bash run_cnrnn.sh
```

To resume from a specific step (e.g., skip Steps 1-2 if already completed):

```bash
bash run_cnrnn.sh --start 3
```

To run only one step:

```bash
bash run_cnrnn.sh --start 4 --end 4
```

## Output

All outputs go to `output/` (configurable via `OUTPUT_DIR` in `run_cnrnn.sh`):

```
output/
├── gc_mapp/
│   └── Annotated_Targets_GRCh38.txt     # GC + mappability per target region
├── read_counts/
│   ├── ReadCount_Matrix.txt              # Regions x samples read count matrix
│   └── ReadCount_Matrix.rds
├── preprocessing/
│   ├── smoothed_L2R.h5                   # Smoothed log2 ratios (HDF5)
│   ├── smoothed_L2R.rds
│   └── CORRseq_results/
│       └── chr*.alpha.*.csv              # Per-chromosome segmentation calls
└── predictions/
    └── Final_CNV_Predictions.csv         # Final output
```

### Final output columns (`Final_CNV_Predictions.csv`)

| Column | Description |
|--------|-------------|
| `sample_idx` | 1-based sample index |
| `sample_id` | Sample identifier |
| `chrom` | Chromosome |
| `start` | CNV start position (bp) |
| `end` | CNV end position (bp) |
| `length_kb` | CNV length in kilobases |
| `num_markers` | Number of exonic target regions spanning the call |
| `type` | `del` (deletion) or `dup` (duplication) |
| `pred_prob` | LSTM prediction probability (0-1) |
| `final_label` | 1 = confirmed CNV, 0 = rejected |
| `prediction_status` | `RNN_Confirmed`, `RNN_Rejected`, or `Kept_Size_TooLarge` |

## Example: Running on an HPC Cluster (SLURM)

`run_cnrnn.sh` is scheduler-agnostic. Wrap it in a SLURM submission script:

```bash
#!/bin/bash
#SBATCH --job-name=cnrnn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=cnrnn_%j.out

module load conda
conda activate cnrnn_env

cd /path/to/CN-RNN
bash run_cnrnn.sh
```

**Resource guidelines:**
- Step 2 (featureCounts) benefits from multiple threads (set `THREADS` in config)
- Step 3 (segmentation) uses parallel R processes (set `R_CORES` in config)
- Step 4 (LSTM) uses GPU if available; falls back to CPU
- Memory: ~2 GB per sample in the read count matrix; 64 GB is sufficient for ~100 samples

## Running Individual Steps

Each script can be called independently:

```bash
# Step 1: GC and mappability annotation
bash scripts/02-calculate_gc_mapp.sh \
    resources/hglft_genome_GRCh38.bed \
    /path/to/reference.fa \
    /path/to/mappability.bw \
    output/gc_mapp/Annotated_Targets_GRCh38.txt

# Step 2: Read counting
Rscript scripts/03-bam_to_readcount.R \
    --sample_matrix sample_matrix.txt \
    --bed_file resources/hglft_genome_GRCh38.bed \
    --output_dir output/read_counts \
    --threads 16

# Step 3: QC, normalization, smoothing, segmentation
Rscript scripts/04-qc_norm_smooth_segment.R \
    --rc_file output/read_counts/ReadCount_Matrix.txt \
    --gc_mapp_file output/gc_mapp/Annotated_Targets_GRCh38.txt \
    --sample_matrix sample_matrix.txt \
    --output_dir output/preprocessing \
    --n_cores 12 \
    --alpha "0.05,0.01"

# Step 4: LSTM prediction
python scripts/05-cnrnn_predict.py \
    --h5_file output/preprocessing/smoothed_L2R.h5 \
    --calls_dir output/preprocessing/CORRseq_results \
    --del_model models/model_del_20250726.keras \
    --dup_model models/model_dup_20250726.keras \
    --output_dir output/predictions
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `THREADS` | 16 | CPU threads for featureCounts (Step 2) |
| `R_CORES` | 12 | Parallel cores for CORRseq segmentation (Step 3) |
| `ALPHA` | 0.05,0.01 | Significance levels for segmentation |
| `MIN_MARKERS` | 5 | Discard CNV calls with fewer markers |
| `MAX_MARKERS` | 120 | Auto-accept CNV calls exceeding this (skip RNN) |

## Citation

*(Manuscript in preparation)*

## License

*(To be added)*
