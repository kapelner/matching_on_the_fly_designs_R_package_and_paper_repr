#!/usr/bin/env bash

# This script sorts comprehensive test results by error_message and opens them in LibreOffice.
# Note: Automatic row freezing is not supported via CSV command line,
# so please press 'Alt+W, F' or 'View -> Freeze Cells -> Freeze First Row' in Calc.

FILES=$(ls package_tests/comprehensive_tests_results_nc_*.csv 2>/dev/null)

if [ -z "$FILES" ]; then
    echo "No comprehensive test result CSV files found in package_tests/"
    exit 1
fi

echo "Sorting CSV files in-place by error_message (errors first)..."

for FILE in $FILES; do
    # Use R to sort by error_message (NA values last) in-place
    Rscript -e "
        suppressPackageStartupMessages(library(data.table))
        dt = fread('$FILE')
        if ('error_message' %in% names(dt)) {
            # Sort with errors first (non-NA first)
            setorder(dt, -error_message, na.last = TRUE)
            fwrite(dt, '$FILE')
        }
    "
done

echo "Opening files in LibreOffice Calc..."
libreoffice --calc --nolockcheck -o package_tests/comprehensive_tests_results_nc_*.csv &
