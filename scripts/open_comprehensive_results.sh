#!/usr/bin/env bash

# This script sorts comprehensive test results by error_message and opens them in LibreOffice.
# Note: Automatic row freezing is not supported via CSV command line,
# so please press 'Alt+W, F' or 'View -> Freeze Cells -> Freeze First Row' in Calc.

FILES=$(ls package_tests/comprehensive_tests_results_nc_*.csv 2>/dev/null)

if [ -z "$FILES" ]; then
    echo "No comprehensive test result CSV files found in package_tests/"
    exit 1
fi

TEMP_DIR=$(mktemp -d)
echo "Preparing sorted files in $TEMP_DIR..."

for FILE in $FILES; do
    BASENAME=$(basename "$FILE")
    SORTED_FILE="$TEMP_DIR/$BASENAME"
    
    # Use R to sort by error_message (NA values last)
    # We use a simple R command to maintain CSV integrity
    Rscript -e "
        library(data.table)
        dt = fread('$FILE')
        if ('error_message' %in% names(dt)) {
            # Sort with errors first (non-NA first)
            setorder(dt, -error_message, na.last = TRUE)
        }
        fwrite(dt, '$SORTED_FILE')
    "
done

echo "Opening sorted files in LibreOffice Calc..."
libreoffice --calc "$TEMP_DIR"/*.csv &

# Optional: cleanup temp dir after a delay or on exit
# (but we want the files to stay until LibreOffice reads them)
(sleep 60 && rm -rf "$TEMP_DIR") &
