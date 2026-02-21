#!/bin/bash
set -e

# ======================================================
# Usage check
# ======================================================
if [ $# -lt 4 ]; then
  echo ""
  echo "Usage:"
  echo "  $0 <PROJECT_DIR> <HYCHPO_MAX> <ABS_HYCH_MAX> <KNOWN_MUTATION_TSV> [BRANCH]"
  echo ""
  echo "Arguments:"
  echo "  PROJECT_DIR          Root directory of the project"
  echo "  HYCHPO_MAX           Max threshold for Hy+Ch-Po (grid search 1~X)"
  echo "  ABS_HYCH_MAX         Max threshold for |Hy-Ch| (grid search 1~Y)"
  echo "  KNOWN_MUTATION_TSV   Ground truth mutation site TSV"
  echo "  BRANCH               Optional: psiblast | smp | all (default: all)"
  echo ""
  echo "Example:"
  echo "  $0 ~/PSSM 10 6 data/known_mutation_sites.tsv all"
  echo ""
  exit 1
fi

# ======================================================
# Input arguments
# ======================================================
PROJECT_DIR=$1
HYCHPO_MAX=$2
ABS_HYCH_MAX=$3
KNOWN_MUTATION_TSV=$4
BRANCH=${5:-all}

# ======================================================
# Validate branch
# ======================================================
if [[ "$BRANCH" != "psiblast" && "$BRANCH" != "smp" && "$BRANCH" != "all" ]]; then
  echo "Error: Invalid branch '$BRANCH'"
  echo "Allowed: psiblast | smp | all"
  exit 1
fi

# ======================================================
# Validate known mutation file
# ======================================================
if [ ! -f "$KNOWN_MUTATION_TSV" ]; then
  echo "Error: known mutation TSV not found:"
  echo "  $KNOWN_MUTATION_TSV"
  echo ""
  echo "The TSV file must contain columns:"
  echo "  ID"
  echo "  Known Mutation Sites"
  echo ""
  exit 1
fi

# ======================================================
# Conda environment
# ======================================================
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate pssm

# ======================================================
# Branch runner
# ======================================================
run_branch_screen() {
  local branch_name=$1

  local input_dir="$PROJECT_DIR/results/pssm/filtered/$branch_name"

  echo ""
  echo "======================================================"
  echo " Running Mutation Site Screening (Grid Search)"
  echo " Branch: $branch_name"
  echo " Input : $input_dir"
  echo "======================================================"
  echo ""

  if [ ! -d "$input_dir" ]; then
    echo "Error: Filtered directory not found:"
    echo "  $input_dir"
    echo ""
    echo "Please run scripts/conservation_filter.sh first."
    exit 1
  fi

  python "$PROJECT_DIR/src/main.py" \
    --stage mutation_site_screen \
    --branch "$branch_name" \
    --hychpo_max "$HYCHPO_MAX" \
    --abs_hych_max "$ABS_HYCH_MAX" \
    --known_mutation_sites_tsv "$KNOWN_MUTATION_TSV"
}

# ======================================================
# Run requested branch(es)
# ======================================================
if [ "$BRANCH" == "all" ]; then
  run_branch_screen "psiblast"
  run_branch_screen "smp"
else
  run_branch_screen "$BRANCH"
fi

# ======================================================
# Deactivate
# ======================================================
conda deactivate

echo ""
echo "✅ Mutation site screening finished."
echo ""