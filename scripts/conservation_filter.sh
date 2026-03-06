#!/bin/bash
set -e

# ======================================================
# Usage check
# ======================================================
if [ $# -lt 5 ]; then
  echo ""
  echo "Usage:"
  echo "  $0 <PROJECT_DIR> <CONS_MIN> <CONS_MAX> <EC_MIN> <EC_MAX> [BRANCH]"
  echo ""
  echo "Arguments:"
  echo "  PROJECT_DIR   Root directory of the project"
  echo "  CONS_MIN      Lower bound for Scorecons conservation"
  echo "  CONS_MAX      Upper bound for Scorecons conservation"
  echo "  EC_MIN        Lower bound for evolutionary conservation"
  echo "  EC_MAX        Upper bound for evolutionary conservation"
  echo "  BRANCH        Optional: psiblast | smp | all (default: all)"
  echo ""
  echo "Example:"
  echo "  $0 ~/PSSM 0.15 1 -1.5 10 all"
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
BRANCH=${6:-all}

# ======================================================
# Validate branch
# ======================================================
if [[ "$BRANCH" != "psiblast" && "$BRANCH" != "smp" && "$BRANCH" != "all" ]]; then
  echo "Error: Invalid branch '$BRANCH'"
  echo "Allowed: psiblast | smp | all"
  exit 1
fi

# ======================================================
# Conda environment
# ======================================================
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate pssm

# ======================================================
# Branch loop
# ======================================================
run_branch_filter() {
  local branch_name=$1

  local integrated_dir="$PROJECT_DIR/results/pssm/integrated/$branch_name"
  local filtered_dir="$PROJECT_DIR/results/pssm/filtered/$branch_name"

  echo ""
  echo "======================================================"
  echo " Running Conservation Filtering Stage"
  echo " Branch: $branch_name"
  echo " Input : $integrated_dir"
  echo " Output: $filtered_dir"
  echo "======================================================"
  echo ""

  if [ ! -d "$integrated_dir" ]; then
    echo "Error: Integrated directory not found:"
    echo "  $integrated_dir"
    echo ""
    echo "Please run scripts/conservation.sh first."
    exit 1
  fi

  python "$PROJECT_DIR/src/main.py" \
    --stage conservation_filter \
    --pssm_integrated_dir "$integrated_dir" \
    --pssm_filtered_output_dir "$filtered_dir" \
    --ec_min "$EC_MIN" \
    --ec_max "$EC_MAX" \
    --cons_min "$CONS_MIN" \
    --cons_max "$CONS_MAX"
}

# ======================================================
# Run requested branch(es)
# ======================================================
if [ "$BRANCH" == "all" ]; then
  run_branch_filter "psiblast"
  run_branch_filter "smp"
else
  run_branch_filter "$BRANCH"
fi

# ======================================================
# Deactivate
# ======================================================
conda deactivate

echo ""
echo "✅ Conservation filtering finished."
echo ""
