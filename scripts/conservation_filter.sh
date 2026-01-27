#!/bin/bash
set -e

# ======================================================
# Usage check
# ======================================================
if [ $# -lt 5 ]; then
  echo ""
  echo "Usage:"
  echo "  $0 <PROJECT_DIR> <EC_MIN> <EC_MAX> <CONS_MIN> <CONS_MAX>"
  echo ""
  echo "Example:"
  echo "  $0 ~/PSSM 0.5 2.5 0.3 1.0"
  echo ""
  exit 1
fi

# ======================================================
# Input arguments
# ======================================================
PROJECT_DIR=$1
CONS_MIN=$2
CONS_MAX=$3
EC_MIN=$4
EC_MAX=$5

# ======================================================
# Path setup
# ======================================================
EVOLUTIONARY_CONSERVATION_DIR="$PROJECT_DIR/results/evolutionary_conservation"
RECONSTRUCT_DIR="$EVOLUTIONARY_CONSERVATION_DIR/reconstruct"

if [ ! -d "$RECONSTRUCT_DIR" ]; then
  echo "Error: reconstruct directory not found:"
  echo "  $RECONSTRUCT_DIR"
  echo ""
  echo "Please run scripts/conservation.sh first."
  exit 1
fi

# ======================================================
# Conda environment
# ======================================================
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate pssm

# ======================================================
# Run conservation filtering stage
# ======================================================
python "$PROJECT_DIR/src/main.py" \
  --stage conservation_filter \
  --conservation_reconstruct_dir "$RECONSTRUCT_DIR" \
  --ec_min "$EC_MIN" \
  --ec_max "$EC_MAX" \
  --cons_min "$CONS_MIN" \
  --cons_max "$CONS_MAX"

# ======================================================
# Deactivate
# ======================================================
conda deactivate
