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
# Path setup (NEW standardized structure)
# ======================================================
CDSEARCH_DIR="$PROJECT_DIR/results/cdsearch"

# ---------------- Branch A (PSI-BLAST) ----------------
PSSM_PROFILES_DIR="$PROJECT_DIR/results/pssm/profiles/psiblast"
PSSM_PROFILES_PSSM_DIR="$PSSM_PROFILES_DIR/pssm_profiles"

PSSM_MATRICES_DIR="$PROJECT_DIR/results/pssm/matrices/psiblast"
PSSM_RECONSTRUCT_DIR="$PROJECT_DIR/results/pssm/reconstruct/psiblast"

# ---------------- Branch B (SMP Parse) ----------------
SMP_MATRICES_DIR="$PROJECT_DIR/results/pssm/matrices/smp"
SMP_RECONSTRUCT_DIR="$PROJECT_DIR/results/pssm/reconstruct/smp"

# ---------------- CDD Database ----------------
CDD_DB="$PROJECT_DIR/blastdb/cdd/Cdd"
CDD_SMP_ROOT="$PROJECT_DIR/blastdb/cdd"

# ======================================================
# Conda environment setup
# ======================================================
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate pssm

# ======================================================
# Branch A - Stage 1: domain PSI-BLAST (generate .pssm profiles)
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage domain_psiblast \
  --domain_fasta_dir "$CDSEARCH_DIR/domains_fasta/hseq_with_gap" \
  --psiblast_output_dir "$PSSM_PROFILES_DIR" \
  --psiblast_blast_db "$CDD_DB" \
  --psiblast_threads 8

# ======================================================
# Branch A - Stage 2: PSSM matrix extraction (profiles -> matrices)
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage pssm_features \
  --pssm_profiles_dir "$PSSM_PROFILES_PSSM_DIR" \
  --pssm_matrix_output_dir "$PSSM_MATRICES_DIR"

# ======================================================
# Branch A - Stage 3: PSSM reconstruction (matrices -> full-length)
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage pssm_reconstruct \
  --pssm_matrix_dir "$PSSM_MATRICES_DIR" \
  --pssm_reconstruct_output_dir "$PSSM_RECONSTRUCT_DIR" \
  --pssm_fasta_path "$INPUT_FASTA" \
  --pssm_cdsearch_table "$CDSEARCH_DIR/cdsearch_top_hits_detailed.tsv"

# ======================================================
# Branch B - Stage 1: SMP parsing (CDD .smp -> matrices)
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage smp_parse \
  --smp_cdsearch_table "$CDSEARCH_DIR/cdsearch_top_hits_detailed.tsv" \
  --smp_root_dir "$CDD_SMP_ROOT" \
  --smp_matrix_output_dir "$SMP_MATRICES_DIR"

# ======================================================
# Branch B - Stage 2: PSSM reconstruction (SMP matrices -> full-length)
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage pssm_reconstruct \
  --pssm_matrix_dir "$SMP_MATRICES_DIR" \
  --pssm_reconstruct_output_dir "$SMP_RECONSTRUCT_DIR" \
  --pssm_fasta_path "$INPUT_FASTA" \
  --pssm_cdsearch_table "$CDSEARCH_DIR/cdsearch_top_hits_detailed.tsv"

# ======================================================
# Deactivate conda environment
# ======================================================
conda deactivate
