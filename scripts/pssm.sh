#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ======================================================
# Check input arguments
# ======================================================
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR> [CUSTOM_FASTA]"
  exit 1
fi

PROJECT_DIR=$1

# Default FASTA file (if user doesn't provide CUSTOM_FASTA)
DEFAULT_FASTA="/home/user_312554028/PSSM/test.fasta"

# CUSTOM_FASTA = 第二個輸入參數（可選）
CUSTOM_FASTA=${2:-""}

# INPUT_FASTA = 若沒給第二參數 → 使用 default
INPUT_FASTA="${CUSTOM_FASTA:-$DEFAULT_FASTA}"

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
# Stage 1: domain_psiblast
# ======================================================
python "$PROJECT_DIR/src/main.py" \
    --stage domain_psiblast \
    --domain_fasta_dir "$OUTPUT_DIR1/domains_fasta/hseq_with_gap" \
    --psiblast_output_dir "$OUTPUT_DIR2" \
    --psiblast_blast_db "$CDD_DB/Cdd" \
    --psiblast_threads 8

# ======================================================
# Stage 2: pssm_features
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage pssm_features \
  --pssm_profiles_dir "$OUTPUT_DIR2/pssm_profiles" \
  --pssm_matrix_output_dir "$OUTPUT_DIR2"

# ======================================================
# Stage 3: pssm_reconstruct
# Use INPUT_FASTA instead of hardcoded FASTA
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage pssm_reconstruct \
  --pssm_fasta_path "$INPUT_FASTA" \
  --pssm_cdsearch_table "$OUTPUT_DIR1/cdsearch_top_hits_detailed.tsv" \
  --pssm_matrix_dir "$OUTPUT_DIR2/pssm_matrices" \
  --pssm_reconstruct_output "$OUTPUT_DIR3"

# ======================================================
# Deactivate environment
# ======================================================
conda deactivate
