#!/bin/bash

# ============================================================================
# CBP2025 Report Generation & Organization Script
# ============================================================================
# This script runs all analysis scripts and organizes outputs into
# well-structured directories for easy navigation and comparison.
# ============================================================================

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 <branch_predictor_name> [path/to/csv]"
    echo "Example: $0 baseline"
    echo "Example: $0 tage-sc-l-alberto-ros results/tage-sc-l-alberto-ros-results.csv"
    exit 2
fi

PREDICTOR_NAME="$1"

# Make a filesystem-safe name for folder + filename prefix
SAFE_NAME=$(echo "$PREDICTOR_NAME" | tr ' /' '__' | tr -cd '[:alnum:]_.-')

if [ -z "$SAFE_NAME" ]; then
    echo "Error: predictor name becomes empty after sanitizing."
    exit 2
fi

# Store project root directory
PROJECT_ROOT="$(pwd)"

# Activate virtual environment if it exists
if [ -d "${PROJECT_ROOT}/.venv" ]; then
    source "${PROJECT_ROOT}/.venv/bin/activate"
fi

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Base directories
RESULTS_CSV="results.csv"  # analysis scripts expect this exact filename
REPORTS_ROOT="${PROJECT_ROOT}/reports/${SAFE_NAME}"
REPORTS_DIR="${REPORTS_ROOT}"
SCRIPT_DIR="scripts/analysis/report_generators"

# Tell the Python scripts where to write outputs, and how to prefix graph filenames
export CBP_REPORTS_DIR="${REPORTS_ROOT}"
export CBP_GRAPH_PREFIX="${SAFE_NAME}__"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  CBP2025 Analysis Pipeline${NC}"
echo -e "${BLUE}  Predictor: ${PREDICTOR_NAME}${NC}"
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

# Resolve which CSV to analyze.
# Priority:
#  1) user-provided path (arg2)
#  2) results/<predictor>.csv or common variants
#  3) legacy results.csv locations
CSV_PATH=""

if [ $# -ge 2 ]; then
    if [ -f "$2" ]; then
        CSV_PATH="$2"
    else
        echo -e "${RED}✗${NC} Provided CSV path not found: $2"
        exit 1
    fi
else
    for candidate in \
        "results/${PREDICTOR_NAME}.csv" \
        "results/${SAFE_NAME}.csv" \
        "results/${PREDICTOR_NAME}-results.csv" \
        "results/${SAFE_NAME}-results.csv" \
        "results/${PREDICTOR_NAME}_results.csv" \
        "results/${SAFE_NAME}_results.csv" \
        "${PREDICTOR_NAME}.csv" \
        "${SAFE_NAME}.csv" \
        "${RESULTS_CSV}" \
        "results/${RESULTS_CSV}" \
        "${SCRIPT_DIR}/${RESULTS_CSV}"; do
        if [ -f "$candidate" ]; then
            CSV_PATH="$candidate"
            break
        fi
    done
fi

if [ -z "$CSV_PATH" ]; then
    echo -e "${RED}✗${NC} No suitable CSV found for predictor '${PREDICTOR_NAME}'."
    echo -e "${YELLOW}Looked for:${NC}"
    echo "  - results/${PREDICTOR_NAME}.csv (and common variants)"
    echo "  - ./results.csv or results/results.csv"
    echo -e "${YELLOW}Hint:${NC} Run traces first or pass an explicit CSV path."
    exit 1
fi

echo -e "   Using: ${CSV_PATH}"

# Copy results.csv to script directory if not already there
cp "${CSV_PATH}" "${SCRIPT_DIR}/${RESULTS_CSV}"

cd ${SCRIPT_DIR}
python3 baselineDifficulty.py > ${REPORTS_DIR}/06_difficulty_analysis/intrinsic_hardness/baseline_difficulty_report.txt 2>&1
cd ${PROJECT_ROOT}
echo -e "${GREEN}✓${NC} Baseline difficulty analysis complete"
echo -e "   Output: ${REPORTS_DIR}/06_difficulty_analysis/intrinsic_hardness/baseline_difficulty_report.txt\n"

# ============================================================================
# Step 3: Analyze phase behavior (MPKI vs 50PercMPKI)
# ============================================================================
echo -e "${YELLOW}[3/6]${NC} Analyzing phase behavior patterns..."

cd ${SCRIPT_DIR}
python3 compareMPKI_vs50PercMPKI.py > ${REPORTS_DIR}/03_phase_behavior/phase_behavior_report.txt 2>&1
cd ${PROJECT_ROOT}

echo -e "${GREEN}✓${NC} Phase behavior analysis complete"
echo -e "   Output: ${REPORTS_DIR}/03_phase_behavior/\n"

# ============================================================================
# Step 3.5: Generate MPKI comparison graphs
# ============================================================================
echo -e "${YELLOW}[3.5/6]${NC} Generating MPKI comparison visualizations..."

cd ${SCRIPT_DIR}
python3 graphComparator.py > ${REPORTS_DIR}/01_misprediction_analysis/mpki_comparison_report.txt 2>&1
cd ${PROJECT_ROOT}

echo -e "${GREEN}✓${NC} MPKI comparison graphs complete"
echo -e "   Output: ${REPORTS_DIR}/01_misprediction_analysis/\n"

# ============================================================================
# Step 4: Compare benchmarks with similar MPKI but different behavior
# ============================================================================
echo -e "${YELLOW}[4/6]${NC} Finding benchmarks with divergent behavior..."

cd ${SCRIPT_DIR}
python3 compareMetrics.py > ${REPORTS_DIR}/05_comparative_analysis/divergent_behavior_report.txt 2>&1
cd ${PROJECT_ROOT}

echo -e "${GREEN}✓${NC} Comparative analysis complete"
echo -e "   Output: ${REPORTS_DIR}/05_comparative_analysis/divergent_behavior_report.txt\n"

# ============================================================================
# Step 5: Generate comprehensive performance graphs
# ============================================================================
echo -e "${YELLOW}[5/6]${NC} Generating performance metric visualizations..."

cd ${SCRIPT_DIR}
python3 generate_graphs.py 2>&1 | tee ${REPORTS_DIR}/04_performance_metrics/graph_generation.log
cd ${PROJECT_ROOT}

echo -e "${GREEN}✓${NC} Graph generation complete"
echo -e "   Output: ${REPORTS_DIR}/04_performance_metrics/\n"

# ============================================================================
# Step 6: Analyze branch density across benchmarks
# ============================================================================
echo -e "${YELLOW}[6/6]${NC} Analyzing branch density patterns..."

cd ${SCRIPT_DIR}
python3 rangeOfBranchesAcrossBenchmarks.py > ${REPORTS_DIR}/02_branch_density/branch_density_report.txt 2>&1
cd ${PROJECT_ROOT}

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
