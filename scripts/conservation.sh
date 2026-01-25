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
CDSEARCH_DIR="$PROJECT_DIR/results/cdsearch_results"
PSSM_RECONSTRUCT_DIR="$PROJECT_DIR/results/domain_psiblast/pssm_reconstruct"
CONSERVATION_DIR="$PROJECT_DIR/results/conservation"

# ======================================================
# Check conservation input directory
# ======================================================
if [ ! -d "$CONSERVATION_DIR" ]; then
  echo "Error: Conservation directory not found:"
  echo "  $CONSERVATION_DIR"
  exit 1
fi

# ======================================================
# Conda environment setup
# ======================================================
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate pssm

# ======================================================
# Conservation reconstruction stage
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage conservation_reconstruct \
  --conservation_fasta_path "$INPUT_FASTA" \
  --conservation_cdsearch_table "$CDSEARCH_DIR/cdsearch_top_hits_detailed.tsv" \
  --pssm_reconstruct_dir "$PSSM_RECONSTRUCT_DIR" \
  --conservation_dir "$CONSERVATION_DIR"

# ======================================================
# Deactivate conda environment
# ======================================================
conda deactivate
