#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Check input arguments
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR>"
  exit 1
fi

# Set paths
PROJECT_DIR=$1
OUTPUT_DIR1="$PROJECT_DIR/results/cdsearch_results"
OUTPUT_DIR2="$PROJECT_DIR/results/domain_psiblast"
OUTPUT_DIR3="$PROJECT_DIR/results/pssm_reconstruct"
CDD_DB="$PROJECT_DIR/blastdb/cdd"
DEFAULT_FASTA="/home/user_312554028/PSSM/test.fasta"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
conda activate pssm

python "$PROJECT_DIR/src/main.py" \
    --stage domain_psiblast \
    --domain_fasta_dir "$OUTPUT_DIR1/domains_fasta/hseq_with_gap" \
    --psiblast_output_dir "$OUTPUT_DIR2" \
    --psiblast_blast_db "$CDD_DB/Cdd" \
    --psiblast_threads 8

python "$PROJECT_DIR/src/main.py" \
  --stage pssm_features \
  --pssm_profiles_dir "$OUTPUT_DIR2/pssm_profiles" \
  --pssm_matrix_output_dir "$OUTPUT_DIR2"

python "$PROJECT_DIR/src/main.py" \
  --stage pssm_reconstruct \
  --pssm_fasta_path "$DEFAULT_FASTA" \
  --pssm_cdsearch_table "$OUTPUT_DIR1/cdsearch_top_hits_detailed.tsv" \
  --pssm_matrix_dir "$OUTPUT_DIR2/pssm_matrices" \
  --pssm_reconstruct_output "$OUTPUT_DIR3"

# Deactivate the Conda environment
conda deactivate