#!/bin/bash

# ============================================================================
# CBP2025 Report Generation & Organization Script
# ============================================================================
# This script runs all analysis scripts and organizes outputs into 
# well-structured directories for easy navigation and comparison.
# ============================================================================

set -e  # Exit on error

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Base directories
RESULTS_CSV="results.csv"  # CSV can be in root or results/ directory
REPORTS_DIR="Reports"
SCRIPT_DIR="ReportGenerators"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  CBP2025 Analysis Pipeline${NC}"
echo -e "${BLUE}========================================${NC}\n"

# ============================================================================
# Step 1: Create organized directory structure
# ============================================================================
echo -e "${YELLOW}[1/6]${NC} Creating directory structure..."

# Main report categories
mkdir -p ${REPORTS_DIR}/01_misprediction_analysis
mkdir -p ${REPORTS_DIR}/02_branch_density
mkdir -p ${REPORTS_DIR}/03_phase_behavior
mkdir -p ${REPORTS_DIR}/04_performance_metrics
mkdir -p ${REPORTS_DIR}/05_comparative_analysis
mkdir -p ${REPORTS_DIR}/06_difficulty_analysis

# Subdirectories for each category
mkdir -p ${REPORTS_DIR}/01_misprediction_analysis/{by_category,overall,graphs}
mkdir -p ${REPORTS_DIR}/02_branch_density/{statistics,distributions,graphs}
mkdir -p ${REPORTS_DIR}/03_phase_behavior/{consistent,variable,graphs}
mkdir -p ${REPORTS_DIR}/04_performance_metrics/{ipc,cycles,graphs}
mkdir -p ${REPORTS_DIR}/05_comparative_analysis/{similar_mpki,divergent_behavior}
mkdir -p ${REPORTS_DIR}/06_difficulty_analysis/{intrinsic_hardness,rankings}

echo -e "${GREEN}✓${NC} Directory structure created\n"

# ============================================================================
# Step 2: Run baseline difficulty analysis
# ============================================================================
echo -e "${YELLOW}[2/6]${NC} Analyzing baseline difficulty (intrinsic hardness)..."

# Check for results.csv in multiple locations
if [ -f "${RESULTS_CSV}" ]; then
    CSV_PATH="${RESULTS_CSV}"
elif [ -f "results/${RESULTS_CSV}" ]; then
    CSV_PATH="results/${RESULTS_CSV}"
elif [ -f "${SCRIPT_DIR}/${RESULTS_CSV}" ]; then
    CSV_PATH="${SCRIPT_DIR}/${RESULTS_CSV}"
else
    echo -e "${RED}✗${NC} results.csv not found in root, results/, or ReportGenerators/."
    echo -e "${YELLOW}Hint:${NC} Run traces first: python scripts/trace_exec_training_list.py --trace_dir traces/ --results_dir results/\n"
    exit 1
fi

echo -e "   Using: ${CSV_PATH}"

# Copy results.csv to ReportGenerators if not already there
if [ ! -f "${SCRIPT_DIR}/${RESULTS_CSV}" ]; then
    cp "${CSV_PATH}" "${SCRIPT_DIR}/${RESULTS_CSV}"
fi

cd ${SCRIPT_DIR}
python baselineDifficulty.py > ../${REPORTS_DIR}/06_difficulty_analysis/intrinsic_hardness/baseline_difficulty_report.txt 2>&1
cd ..
echo -e "${GREEN}✓${NC} Baseline difficulty analysis complete"
echo -e "   Output: ${REPORTS_DIR}/06_difficulty_analysis/intrinsic_hardness/baseline_difficulty_report.txt\n"

# ============================================================================
# Step 3: Analyze phase behavior (MPKI vs 50PercMPKI)
# ============================================================================
echo -e "${YELLOW}[3/6]${NC} Analyzing phase behavior patterns..."

cd ${SCRIPT_DIR}
python compareMPKI_vs50PercMPKI.py > ../${REPORTS_DIR}/03_phase_behavior/phase_behavior_report.txt 2>&1
cd ..

echo -e "${GREEN}✓${NC} Phase behavior analysis complete"
echo -e "   Output: ${REPORTS_DIR}/03_phase_behavior/\n"

# ============================================================================
# Step 3.5: Generate MPKI comparison graphs
# ============================================================================
echo -e "${YELLOW}[3.5/6]${NC} Generating MPKI comparison visualizations..."

cd ${SCRIPT_DIR}
python graphComparator.py > ../${REPORTS_DIR}/01_misprediction_analysis/mpki_comparison_report.txt 2>&1
cd ..

echo -e "${GREEN}✓${NC} MPKI comparison graphs complete"
echo -e "   Output: ${REPORTS_DIR}/01_misprediction_analysis/\n"

# ============================================================================
# Step 4: Compare benchmarks with similar MPKI but different behavior
# ============================================================================
echo -e "${YELLOW}[4/6]${NC} Finding benchmarks with divergent behavior..."

cd ${SCRIPT_DIR}
python compareMetrics.py > ../${REPORTS_DIR}/05_comparative_analysis/divergent_behavior_report.txt 2>&1
cd ..

echo -e "${GREEN}✓${NC} Comparative analysis complete"
echo -e "   Output: ${REPORTS_DIR}/05_comparative_analysis/divergent_behavior_report.txt\n"

# ============================================================================
# Step 5: Generate comprehensive performance graphs
# ============================================================================
echo -e "${YELLOW}[5/6]${NC} Generating performance metric visualizations..."

cd ${SCRIPT_DIR}
python generate_graphs.py 2>&1 | tee ../${REPORTS_DIR}/04_performance_metrics/graph_generation.log
cd ..

echo -e "${GREEN}✓${NC} Graph generation complete"
echo -e "   Output: ${REPORTS_DIR}/04_performance_metrics/\n"

# ============================================================================
# Step 6: Analyze branch density across benchmarks
# ============================================================================
echo -e "${YELLOW}[6/6]${NC} Analyzing branch density patterns..."

cd ${SCRIPT_DIR}
python rangeOfBranchesAcrossBenchmarks.py > ../${REPORTS_DIR}/02_branch_density/branch_density_report.txt 2>&1
cd ..

# Move branch density graphs
if ls *.png 1> /dev/null 2>&1; then
    mv *branch*.png ${REPORTS_DIR}/02_branch_density/graphs/ 2>/dev/null || true
    mv *density*.png ${REPORTS_DIR}/02_branch_density/graphs/ 2>/dev/null || true
fi

echo -e "${GREEN}✓${NC} Branch density analysis complete"
echo -e "   Output: ${REPORTS_DIR}/02_branch_density/\n"

# ============================================================================
# Final Summary
# ============================================================================
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  Analysis Complete!${NC}"
echo -e "${BLUE}========================================${NC}\n"

echo -e "${GREEN}All reports organized in: ${REPORTS_DIR}/${NC}\n"

echo "Directory Structure:"
echo "│   ├── statistics/         # Branch per cycle stats"
echo "│   ├── distributions/      # Density distributions"
echo "│   └── graphs/             # Branch density plots"
echo "│"
echo "├── 03_phase_behavior/"
echo "│   ├── consistent/         # Benchmarks with steady behavior"
echo "│   ├── variable/           # Benchmarks with phase changes"
echo "│   └── graphs/             # Phase behavior visualizations"
echo "│"
echo "├── 04_performance_metrics/"
echo "│   ├── ipc/                # Instructions per cycle analysis"
echo "│   ├── cycles/             # Execution cycle metrics"
echo "│   └── graphs/             # Performance visualizations"
echo "│"
echo "├── 05_comparative_analysis/"
echo "│   ├── similar_mpki/       # Benchmarks with similar difficulty"
echo "│   └── divergent_behavior/ # Same MPKI, different characteristics"
echo "│"
echo "└── 06_difficulty_analysis/"
echo "    ├── intrinsic_hardness/ # Baseline difficulty rankings"
echo "    └── rankings/           # Hardest/easiest benchmarks"
echo ""

echo -e "${YELLOW}Key Reports:${NC}"
echo "  • Difficulty: ${REPORTS_DIR}/06_difficulty_analysis/intrinsic_hardness/baseline_difficulty_report.txt"
echo "  • Phase Behavior: ${REPORTS_DIR}/03_phase_behavior/phase_behavior_report.txt"
echo "  • Comparisons: ${REPORTS_DIR}/05_comparative_analysis/divergent_behavior_report.txt"
echo "  • Branch Density: ${REPORTS_DIR}/02_branch_density/branch_density_report.txt"
echo ""

echo -e "${GREEN}✓ Pipeline complete!${NC}"
