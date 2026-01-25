#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ======================================================
# Check input arguments
# ======================================================
if [ $# -lt 2 ]; then
  echo "Usage: $0 <PROJECT_DIR> <INPUT_FASTA>"
  echo "Error: An input FASTA file is required."
  exit 1
fi

# ======================================================
# Set input paths
# ======================================================
PROJECT_DIR=$1
INPUT_FASTA=$2

# Verify that the input FASTA file exists
if [ ! -f "$INPUT_FASTA" ]; then
  echo "Error: Input FASTA file not found:"
  echo "  $INPUT_FASTA"
  exit 1
fi

# ======================================================
# Path setup
# ======================================================
OUTPUT_DIR1="$PROJECT_DIR/results/cdsearch_results"
OUTPUT_DIR2="$PROJECT_DIR/results/domain_psiblast"
OUTPUT_DIR3="$PROJECT_DIR/results/pssm_reconstruct"
CDD_DB="$PROJECT_DIR/blastdb/cdd"

# ======================================================
# Conda environment setup
# ======================================================
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate pssm

# ======================================================
# Stage 1: domain PSI-BLAST
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage domain_psiblast \
  --domain_fasta_dir "$OUTPUT_DIR1/domains_fasta/hseq_with_gap" \
  --psiblast_output_dir "$OUTPUT_DIR2" \
  --psiblast_blast_db "$CDD_DB/Cdd" \
  --psiblast_threads 8

# ======================================================
# Stage 2: PSSM feature extraction
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage pssm_features \
  --pssm_root_dir "$OUTPUT_DIR2"

# ======================================================
# Stage 3: PSSM reconstruction
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage pssm_reconstruct \
  --pssm_root_dir "$OUTPUT_DIR2" \
  --pssm_fasta_path "$INPUT_FASTA" \
  --pssm_cdsearch_table "$OUTPUT_DIR1/cdsearch_top_hits_detailed.tsv"

# ======================================================
# Deactivate conda environment
# ======================================================
conda deactivate
