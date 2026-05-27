#!/usr/bin/env bash
# Searches for .cram files matching sample IDs of interest in
# ./00_Max_samples and ./01_Max_samples in Sophia's disk.

OUTPUT_FILE="/Users/nv4/gitClones/08_deconvTable/083_Results/08_deconvTable_sophiaArchivedSamples/DFT1_cram_paths.txt"

# Sample ID DFT1_PREFIXES to search for
# Pattern: everything up to and including the capital T (e.g. 1058T matches 1058T1a, 1058T2, etc.)
DFT1_PREFIXES=(
    # Diploid
    "84T"
    "134T"
    "142T"
    "146T"
    "1011T"
    "1058T"
    "1059T"
    # Tetraploid
    "18T"
    "477T"
    "528T"
)

SEARCH_DIRS=(
    "./00_Max_samples"
    "./01_Max_samples"
)

echo "Searching for .cram files matching sample IDs of interest..."
echo "Search directories: ${SEARCH_DIRS[*]}"
echo ""

# Build the find command with -name patterns for each prefix
FIND_ARGS=()
for i in "${!DFT1_PREFIXES[@]}"; do
    prefix="${DFT1_PREFIXES[$i]}"
    if [ $i -eq 0 ]; then
        FIND_ARGS+=( -name "${prefix}*.cram" )
    else
        FIND_ARGS+=( -o -name "${prefix}*.cram" )
    fi
done

# Run find across both directories
MATCHES=$(find "${SEARCH_DIRS[@]}" \
    -type f \( "${FIND_ARGS[@]}" \) \
    2>/dev/null | sort)

if [ -z "$MATCHES" ]; then
    echo "No matching .cram files found."
    exit 0
fi

# Print and save results
echo "$MATCHES" | tee "$OUTPUT_FILE"

COUNT=$(echo "$MATCHES" | wc -l)
echo ""
echo "Found $COUNT matching .cram file(s)."
echo "Full paths saved to: $OUTPUT_FILE"