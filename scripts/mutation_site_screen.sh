#!/bin/bash
set -e

# ======================================================
# Usage check
# ======================================================
if [ $# -lt 11 ]; then
  echo ""
  echo "Usage:"
  echo "  $0 <PROJECT_DIR> <X_MIN> <X_MAX> <Y_MIN> <Y_MAX> <Z_MIN> <Z_MAX> <SCORE_THR_MAX> <FEATURE_THR_MAX> <KNOWN_MUTATION_TSV> <N_JOBS> [BRANCH]"
  echo ""
  echo "Arguments:"
  echo "  PROJECT_DIR          Root directory of the project"
  echo "  X_MIN                Minimum weight for Hy"
  echo "  X_MAX                Maximum weight for Hy"
  echo "  Y_MIN                Minimum weight for Ch"
  echo "  Y_MAX                Maximum weight for Ch"
  echo "  Z_MIN                Minimum weight for Po"
  echo "  Z_MAX                Maximum weight for Po"
  echo "  SCORE_THR_MAX        Maximum threshold for score filtering"
  echo "  FEATURE_THR_MAX      Maximum threshold for feature filtering"
  echo "  KNOWN_MUTATION_TSV   Ground truth mutation site TSV"
  echo "  N_JOBS               Number of CPU cores"
  echo "  BRANCH               Optional: psiblast | smp | all (default: all)"
  echo ""
  echo "Example:"
  echo "  $0 ~/PSSM -10 10 -10 10 -10 10 10 10 data/known_mutation_sites.tsv 32 all"
  echo ""
  exit 1
fi

# ======================================================
# Input arguments
# ======================================================
PROJECT_DIR=$1
X_MIN=$2
X_MAX=$3
Y_MIN=$4
Y_MAX=$5
Z_MIN=$6
Z_MAX=$7
SCORE_THR_MAX=$8
FEATURE_THR_MAX=$9
KNOWN_MUTATION_TSV=${10}
N_JOBS=${11}
BRANCH=${12:-all}

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
# Validate CPU cores
# ======================================================
if ! [[ "$N_JOBS" =~ ^[0-9]+$ ]]; then
  echo "Error: N_JOBS must be an integer"
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
  echo " Running Mutation Site Screening (Ultra Fast Grid Search)"
  echo " Branch : $branch_name"
  echo " Input  : $input_dir"
  echo " CPUs   : $N_JOBS"
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
    --x_min "$X_MIN" \
    --x_max "$X_MAX" \
    --y_min "$Y_MIN" \
    --y_max "$Y_MAX" \
    --z_min "$Z_MIN" \
    --z_max "$Z_MAX" \
    --score_thr_max "$SCORE_THR_MAX" \
    --feature_thr_max "$FEATURE_THR_MAX" \
    --known_mutation_sites_tsv "$KNOWN_MUTATION_TSV" \
    --n_jobs "$N_JOBS"
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