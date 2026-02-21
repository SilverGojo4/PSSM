#!/bin/bash

set -e

# ======================================================
# Usage
# ======================================================
if [ $# -lt 2 ]; then
  echo "Usage: $0 <PROJECT_DIR> <INPUT_FASTA> [BRANCH]"
  echo ""
  echo "Arguments:"
  echo "  PROJECT_DIR   Root directory of the project"
  echo "  INPUT_FASTA   Input FASTA file (full-length sequences)"
  echo "  BRANCH        Optional: psiblast | smp | all (default: all)"
  echo ""
  exit 1
fi

# ======================================================
# Input args
# ======================================================
PROJECT_DIR=$1
INPUT_FASTA=$2
BRANCH=${3:-all}

# ======================================================
# Validate branch
# ======================================================
if [[ "$BRANCH" != "psiblast" && "$BRANCH" != "smp" && "$BRANCH" != "all" ]]; then
  echo "Error: Invalid branch '$BRANCH'"
  echo "Allowed: psiblast | smp | all"
  exit 1
fi

# ======================================================
# Verify FASTA
# ======================================================
if [ ! -f "$INPUT_FASTA" ]; then
  echo "Error: Input FASTA file not found:"
  echo "  $INPUT_FASTA"
  exit 1
fi

# ======================================================
# Shared paths
# ======================================================
CDSEARCH_TABLE="$PROJECT_DIR/results/cdsearch/cdsearch_top_hits_detailed.tsv"
SCORECONS_DIR="$PROJECT_DIR/results/scorecons"
CONSURF_DIR="$PROJECT_DIR/results/consurf"

# ======================================================
# Check required dirs/files
# ======================================================
if [ ! -f "$CDSEARCH_TABLE" ]; then
  echo "Error: CD-Search table not found:"
  echo "  $CDSEARCH_TABLE"
  exit 1
fi

if [ ! -d "$SCORECONS_DIR" ]; then
  echo "Error: Scorecons directory not found:"
  echo "  $SCORECONS_DIR"
  exit 1
fi

if [ ! -d "$CONSURF_DIR" ]; then
  echo "Error: ConSurf directory not found:"
  echo "  $CONSURF_DIR"
  exit 1
fi

# ======================================================
# Conda environment setup
# ======================================================
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate pssm

# ======================================================
# Branch loop
# ======================================================
run_branch_pipeline() {
  local branch_name=$1
  local integrated_dir="$PROJECT_DIR/results/pssm/integrated/$branch_name"

  echo ""
  echo "======================================================"
  echo " Running Conservation Integration Pipeline"
  echo " Branch: $branch_name"
  echo " Integrated Dir: $integrated_dir"
  echo "======================================================"
  echo ""

  # ------------------------------------------------------
  # Stage 1: Scorecons reconstruction
  # ------------------------------------------------------
  python "$PROJECT_DIR/src/main.py" \
    --stage scorecons_reconstruct \
    --branch "$branch_name" \
    --pssm_fasta_path "$INPUT_FASTA" \
    --pssm_cdsearch_table "$CDSEARCH_TABLE" \
    --scorecons_dir "$SCORECONS_DIR" \
    --pssm_integrated_output_dir "$integrated_dir"

  # ------------------------------------------------------
  # Stage 2: ConSurf integration
  # ------------------------------------------------------
  python "$PROJECT_DIR/src/main.py" \
    --stage consurf_integrate \
    --branch "$branch_name" \
    --pssm_fasta_path "$INPUT_FASTA" \
    --consurf_dir "$CONSURF_DIR" \
    --pssm_integrated_dir "$integrated_dir"
}

# ======================================================
# Run requested branch(es)
# ======================================================
if [ "$BRANCH" == "all" ]; then
  run_branch_pipeline "psiblast"
  run_branch_pipeline "smp"
else
  run_branch_pipeline "$BRANCH"
fi

# ======================================================
# Deactivate conda environment
# ======================================================
conda deactivate

echo ""
echo "✅ Conservation integration pipeline finished."
echo ""
