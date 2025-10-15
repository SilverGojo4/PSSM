#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Check input arguments
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR> [CUSTOM_FASTA]"
  exit 1
fi

# Set paths
PROJECT_DIR=$1
CUSTOM_FASTA=${2:-""}
DEFAULT_FASTA="$PROJECT_DIR/data/processed/proteins_wt.fasta"
INPUT_FASTA="${CUSTOM_FASTA:-$DEFAULT_FASTA}"
OUTPUT_DIR1="$PROJECT_DIR/results/cdsearch_results"
CDD_DB="$PROJECT_DIR/blastdb/cdd"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
conda activate pssm

python "$PROJECT_DIR/src/main.py" \
  --stage cdsearch \
  --cdsearch_input_fasta "$INPUT_FASTA" \
  --cdsearch_output_dir "$OUTPUT_DIR1" \
  --cdsearch_cdd_db "$CDD_DB"

# Deactivate the Conda environment
conda deactivate