#!/bin/bash

# Backwards-compatible wrapper.
# The pipeline logic moved to generateReportForPredictor.sh so it's easier to find.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec "${SCRIPT_DIR}/generateReportForPredictor.sh" "$@"
