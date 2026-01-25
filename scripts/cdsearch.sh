#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ======================
# Check input arguments
# ======================
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR> [CUSTOM_FASTA]"
  exit 1
fi

# ======================
# Set input paths
# ======================
PROJECT_DIR=$1
INPUT_FASTA=$2

# Check whether the input FASTA file exists
if [ ! -f "$INPUT_FASTA" ]; then
  echo "Error: Input FASTA file not found:"
  echo "  $INPUT_FASTA"
  exit 1
fi

# ======================
# Output and database paths
# ======================
OUTPUT_DIR1="$PROJECT_DIR/results/cdsearch_results"
CDD_DB="$PROJECT_DIR/blastdb/cdd"
CDD_DB_TABLE="$PROJECT_DIR/blastdb/cddid_all.tbl"

# ======================
# Conda environment setup
# ======================
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate pssm

# ======================
# Run Python pipeline
# ======================
python "$PROJECT_DIR/src/main.py" \
  --stage cdsearch \
  --cdsearch_input_fasta "$INPUT_FASTA" \
  --cdsearch_output_dir "$OUTPUT_DIR1" \
  --cdsearch_cdd_db "$CDD_DB" \
  --cdsearch_cddid_tbl "$CDD_DB_TABLE"

# ======================
# Deactivate conda environment
# ======================
conda deactivate
