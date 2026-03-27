#!/bin/bash

# Script to run all submissions sequentially
# This MUST run sequentially because each submission overwrites the same source files

set -e  # Exit on error

# Check if running in tmux, if not, start a tmux session
if [ -z "$TMUX" ]; then
    # Check if tmux is installed
    if ! command -v tmux &> /dev/null; then
        echo "WARNING: tmux is not installed. Running without tmux..."
        echo "To install tmux: brew install tmux (macOS) or apt-get install tmux (Linux)"
        echo ""
        sleep 2
        # Continue without tmux
    else
        SESSION_NAME="cbp-submissions-$(date +%Y%m%d-%H%M%S)"
        echo "Starting tmux session: $SESSION_NAME"
        echo "You can detach with Ctrl+B then D, and reattach with: tmux attach -t $SESSION_NAME"
        echo ""
        sleep 2
        tmux new-session -s "$SESSION_NAME" "$0 --in-tmux; echo ''; echo 'Press Enter to close tmux session...'; read"
        exit 0
    fi
fi

# Parse arguments
IN_TMUX=false
if [[ "$1" == "--in-tmux" ]]; then
    IN_TMUX=true
fi

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
    # Remove results directory
    rm -rf ./results/baseline/ 2>/dev/null || true

    # Remove build artifacts
    make clean >/dev/null 2>&1 || true

    # Remove ALL .h and .cc files (they're submission-specific)
    find . -maxdepth 1 -name "*.h" -delete 2>/dev/null || true
    rm -f *.cc 2>/dev/null || true

    # Restore minimal base framework file needed by Makefile
    # (my_cond_branch_predictor.cc can be empty - submissions provide .h)
    cp BaseFramework/my_cond_branch_predictor.cc . 2>/dev/null || touch my_cond_branch_predictor.cc
}

# Function to draw progress bar
draw_progress_bar() {
    local current=$1
    local total=$2
    local width=50
    local percentage=$((100 * current / total))
    local filled=$((width * current / total))
    local empty=$((width - filled))

    printf "\r${BLUE}Overall Progress: [${NC}"
    printf "%${filled}s" | tr ' ' '='
    printf "%${empty}s" | tr ' ' '-'
    printf "${BLUE}] ${GREEN}%d/%d${NC} ${YELLOW}(%d%%)${NC}" "$current" "$total" "$percentage"
}

# Function to process a single submission
process_submission() {
    local submission_dir=$1
    local submission_name=$(basename "$submission_dir")
    local current_num=$2
    local total_num=$3

    echo -e "\n\n${BLUE}========================================${NC}"
    echo -e "${BLUE}  SUBMISSION [$current_num/$total_num]: $submission_name${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo "Processing: $submission_name" >> "$LOG_FILE"

    # Track start time
    submission_start_time=$(date +%s)

    # Skip if not a directory or is .DS_Store
    if [[ ! -d "$submission_dir" ]] || [[ "$submission_name" == ".DS_Store" ]]; then
        echo -e "${YELLOW}Skipping $submission_name (not a valid submission)${NC}"
        return
    fi

    # Step 1: Clean up ALL files from previous submission
    echo -e "${YELLOW}Step 1/4: Cleaning up previous submission...${NC}"
    cleanup
    echo -e "${GREEN}Cleanup complete${NC}"

    # Step 2: Copy all files from submission to main directory
    echo -e "${YELLOW}Step 2/4: Copying files from submission...${NC}"
    cp -r "$submission_dir"/* "$MAIN_DIR/" 2>/dev/null || true
    # Remove .DS_Store if copied
    find "$MAIN_DIR" -maxdepth 1 -name ".DS_Store" -delete 2>/dev/null || true
    echo -e "${GREEN}Files copied${NC}"

    # Step 3: Run the scripts
    echo -e "${YELLOW}Step 3/4: Running ./runScripts.sh...${NC}"
    if ./runScripts.sh 2>&1 | tee -a "$LOG_FILE"; then
        echo -e "${GREEN}runScripts.sh completed successfully${NC}"
    else
        echo -e "${RED}ERROR: runScripts.sh failed for $submission_name${NC}"
        echo "ERROR: runScripts.sh failed for $submission_name" >> "$LOG_FILE"
        cleanup
        return 1
    fi

    # Step 4: Save results
    echo -e "${YELLOW}Step 4/4: Saving results...${NC}"
    if [[ -f "./results/baseline/results.csv" ]]; then
        cp "./results/baseline/results.csv" "$RESULTS_CSV_DIR/${submission_name}.csv"
        echo -e "${GREEN}Results saved to $RESULTS_CSV_DIR/${submission_name}.csv${NC}"
        echo "SUCCESS: Results saved for $submission_name" >> "$LOG_FILE"
    else
        echo -e "${RED}WARNING: results.csv not found for $submission_name${NC}"
        echo "WARNING: results.csv not found for $submission_name" >> "$LOG_FILE"
    fi

    # Calculate and display time taken
    submission_end_time=$(date +%s)
    submission_duration=$((submission_end_time - submission_start_time))
    submission_minutes=$((submission_duration / 60))
    submission_seconds=$((submission_duration % 60))

    echo -e "${GREEN}✓ Completed: $submission_name${NC}"
    echo -e "${BLUE}  Time taken: ${submission_minutes}m ${submission_seconds}s${NC}\n"
    echo "Time taken for $submission_name: ${submission_minutes}m ${submission_seconds}s" >> "$LOG_FILE"
}

# Main execution
if [ "$IN_TMUX" = true ]; then
    echo -e "${GREEN}✓ Running in tmux session${NC}"
    echo -e "${YELLOW}  You can detach with: Ctrl+B then D${NC}"
    echo -e "${YELLOW}  To reattach later: tmux attach${NC}"
    echo ""
fi

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  Running All Submissions Sequentially  ${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "${YELLOW}NOTE: Running sequentially to avoid file conflicts${NC}"
echo ""

# Overall start time
overall_start_time=$(date +%s)

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

        # Draw overall progress bar
        draw_progress_bar "$current" "$total_submissions"

        if process_submission "$submission_dir" "$current" "$total_submissions"; then
            succeeded=$((succeeded + 1))
        else
            failed=$((failed + 1))
        fi

        # Update progress bar after completion
        echo ""  # New line after progress bar
    fi
done

# Clear progress bar line
echo ""

# Calculate overall time
overall_end_time=$(date +%s)
overall_duration=$((overall_end_time - overall_start_time))
overall_hours=$((overall_duration / 3600))
overall_minutes=$(((overall_duration % 3600) / 60))
overall_seconds=$((overall_duration % 60))

# Final summary
echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}         FINAL SUMMARY                  ${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}✓ Succeeded: $succeeded${NC}"
if [[ $failed -gt 0 ]]; then
    echo -e "${RED}✗ Failed: $failed${NC}"
else
    echo -e "${GREEN}✓ Failed: $failed${NC}"
fi
echo -e "${BLUE}  Total: $total_submissions${NC}"
echo ""
echo -e "${BLUE}⏱  Total Time: ${overall_hours}h ${overall_minutes}m ${overall_seconds}s${NC}"
echo ""
echo -e "${GREEN}All results saved in: $RESULTS_CSV_DIR${NC}"
echo -e "${GREEN}Log file: $LOG_FILE${NC}"
echo ""
echo "=== Completed at $(date) ===" >> "$LOG_FILE"
echo "Succeeded: $succeeded, Failed: $failed, Total: $total_submissions" >> "$LOG_FILE"
echo "Total execution time: ${overall_hours}h ${overall_minutes}m ${overall_seconds}s" >> "$LOG_FILE"

# List all generated CSV files
echo -e "${BLUE}Generated CSV files:${NC}"
ls -lh "$RESULTS_CSV_DIR"/*.csv 2>/dev/null || echo "No CSV files generated"

# Tmux session completion message
if [ "$IN_TMUX" = true ]; then
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  All submissions completed!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo -e "${YELLOW}This tmux session will remain open.${NC}"
    echo -e "${YELLOW}To exit tmux: Press Ctrl+B then type :kill-session${NC}"
    echo -e "${YELLOW}Or just close this terminal.${NC}"
fi
