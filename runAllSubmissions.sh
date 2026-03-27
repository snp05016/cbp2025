#!/bin/bash

# Script to run all submissions sequentially
# This MUST run sequentially because each submission overwrites the same source files

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

SUBMISSIONS_DIR="./Submissions"
RESULTS_CSV_DIR="./megascriptruncsvs"
MAIN_DIR=$(pwd)

# Create megascriptruncsvs directory if it doesn't exist
mkdir -p "$RESULTS_CSV_DIR"

# Log file
LOG_FILE="$RESULTS_CSV_DIR/run_log.txt"
echo "=== Starting submission runs at $(date) ===" > "$LOG_FILE"

# Function to clean up after each run
cleanup() {
    echo -e "${YELLOW}Cleaning up...${NC}"
    # Remove results directory
    rm -rf ./results/baseline/
    # Remove build artifacts
    make clean 2>/dev/null || true
    echo -e "${GREEN}Cleanup complete${NC}"
}

# Function to process a single submission
process_submission() {
    local submission_dir=$1
    local submission_name=$(basename "$submission_dir")

    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}Processing: $submission_name${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo "Processing: $submission_name" >> "$LOG_FILE"

    # Skip if not a directory or is .DS_Store
    if [[ ! -d "$submission_dir" ]] || [[ "$submission_name" == ".DS_Store" ]]; then
        echo -e "${YELLOW}Skipping $submission_name (not a valid submission)${NC}"
        return
    fi

    # Step 1: Copy all files from submission to main directory
    echo -e "${YELLOW}Step 1: Copying files from submission...${NC}"
    cp -r "$submission_dir"/* "$MAIN_DIR/" 2>/dev/null || true
    # Remove .DS_Store if copied
    find "$MAIN_DIR" -maxdepth 1 -name ".DS_Store" -delete 2>/dev/null || true
    echo -e "${GREEN}Files copied${NC}"

    # Step 2: Run the scripts
    echo -e "${YELLOW}Step 2: Running ./runScripts.sh...${NC}"
    if ./runScripts.sh 2>&1 | tee -a "$LOG_FILE"; then
        echo -e "${GREEN}runScripts.sh completed successfully${NC}"
    else
        echo -e "${RED}ERROR: runScripts.sh failed for $submission_name${NC}"
        echo "ERROR: runScripts.sh failed for $submission_name" >> "$LOG_FILE"
        cleanup
        return 1
    fi

    # Step 3: Move results.csv to megascriptruncsvs with submission name
    echo -e "${YELLOW}Step 3: Saving results...${NC}"
    if [[ -f "./results/baseline/results.csv" ]]; then
        cp "./results/baseline/results.csv" "$RESULTS_CSV_DIR/${submission_name}.csv"
        echo -e "${GREEN}Results saved to $RESULTS_CSV_DIR/${submission_name}.csv${NC}"
        echo "SUCCESS: Results saved for $submission_name" >> "$LOG_FILE"
    else
        echo -e "${RED}WARNING: results.csv not found for $submission_name${NC}"
        echo "WARNING: results.csv not found for $submission_name" >> "$LOG_FILE"
    fi

    # Step 4: Cleanup
    cleanup

    echo -e "${GREEN}Completed: $submission_name${NC}\n"
}

# Main execution
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  Running All Submissions Sequentially  ${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "${YELLOW}NOTE: Running sequentially to avoid file conflicts${NC}"
echo ""

# Check if Submissions directory exists
if [[ ! -d "$SUBMISSIONS_DIR" ]]; then
    echo -e "${RED}ERROR: Submissions directory not found!${NC}"
    exit 1
fi

# Count total submissions
total_submissions=$(find "$SUBMISSIONS_DIR" -mindepth 1 -maxdepth 1 -type d ! -name ".DS_Store" | wc -l | tr -d ' ')
echo -e "${BLUE}Found $total_submissions submissions to process${NC}\n"

# Counter
current=0
failed=0
succeeded=0

# Process each submission
for submission_dir in "$SUBMISSIONS_DIR"/*; do
    if [[ -d "$submission_dir" ]] && [[ $(basename "$submission_dir") != ".DS_Store" ]]; then
        current=$((current + 1))
        echo -e "${BLUE}[$current/$total_submissions]${NC}"

        if process_submission "$submission_dir"; then
            succeeded=$((succeeded + 1))
        else
            failed=$((failed + 1))
        fi
    fi
done

# Final summary
echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}         FINAL SUMMARY                  ${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}Succeeded: $succeeded${NC}"
echo -e "${RED}Failed: $failed${NC}"
echo -e "${BLUE}Total: $total_submissions${NC}"
echo ""
echo -e "${GREEN}All results saved in: $RESULTS_CSV_DIR${NC}"
echo -e "${GREEN}Log file: $LOG_FILE${NC}"
echo ""
echo "=== Completed at $(date) ===" >> "$LOG_FILE"
echo "Succeeded: $succeeded, Failed: $failed, Total: $total_submissions" >> "$LOG_FILE"

# List all generated CSV files
echo -e "${BLUE}Generated CSV files:${NC}"
ls -lh "$RESULTS_CSV_DIR"/*.csv 2>/dev/null || echo "No CSV files generated"
