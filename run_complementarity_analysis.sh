#!/bin/bash
# Runner script for predictor complementarity analysis

echo "Running Predictor Complementarity Analysis..."
echo "=============================================="
echo

# Make sure we're in the right directory
cd "$(dirname "$0")"

# Run the analysis script
python3 scripts/analysis/predictor_complementarity_analysis.py

echo
echo "Done! Check scripts/analysis/complementarity_analysis/ for results."
