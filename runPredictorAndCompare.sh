#!/bin/bash

# ============================================================================
# CBP2025: Run a predictor report + multi-predictor comparisons
# ============================================================================
# 1) Runs the single-predictor pipeline (generateReportForPredictor.sh)
# 2) Runs scripts/analysis/compare_predictors.py across results/*.csv
# ============================================================================

set -e

if [ $# -lt 1 ]; then
  echo "Usage: $0 <branch_predictor_name> [path/to/csv]"
  echo "Example: $0 baseline"
  echo "Example: $0 tage-sc-l-alberto-ros results/tage-sc-l-alberto-ros-results.csv"
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Prefer repo venv python if it exists
PYTHON_BIN="python3"
if [ -x "${REPO_ROOT}/.venv/bin/python" ]; then
  PYTHON_BIN="${REPO_ROOT}/.venv/bin/python"
fi

# 1) Per-predictor report
"${REPO_ROOT}/generateReportForPredictor.sh" "$@"

# 2) Multi-predictor comparisons (compares all results/*.csv)
echo ""
echo "========================================"
echo " Running multi-predictor comparisons"
echo "========================================"
"${PYTHON_BIN}" "${REPO_ROOT}/scripts/analysis/compare_predictors.py"

echo ""
echo "✓ Done. Multi-predictor outputs: ${REPO_ROOT}/reports/comparison-predictors/"
