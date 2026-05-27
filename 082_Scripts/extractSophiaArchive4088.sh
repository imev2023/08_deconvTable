#!/usr/bin/env bash
# extract_regions.sh
# For each CRAM file of interest, extracts regions of interest,
# indexes them, and saves them organised per region.

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
DATADIR="/Users/nv4/gitClones/08_deconvTable/083_Results/08_deconvTable_sophiaArchivedSamples/DFT1_cram_paths.txt"   # Pass a directory OR use the paths file from find_cram_files.sh
RESDIR="/Users/nv4/gitClones/08_deconvTable/083_Results/08_deconvTable_sophiaArchivedSamples/dft1Germline/" # Output base directory
REF="/Users/nv4/gitClones/08_deconvTable/081_Data/fromFarm/Sarcophilus_harrisii.mSarHar_split.fa"

# # ── Regions of interest ───────────────────────────────────────────────────────
# # Format: "LABEL:CHR:START-END". Of note, these have been aligned to the split reference genome, so regions on CHR1-4 are split by p and q arms!
# # 1p	286703773
# # 1q	429709856
# # 2q	341385470
# # 2p	321366317
# # 3p	192124700
# # 3q	419222568
# # 4	464895054
# # 5	288121652
# # 6	254895979
# # X	83081154
# # Y	130564
# # MT	17117
# REGIONS=(
#     "CHR1_240000000_255000000:1p:240000000-255000000"
#     ## p = 1p:1-286,703,773
#     "CHR1_530000000_600000000:1q:243296227-313296227"

#     ## 2p:1-321,366,317
#     # "CHR2_240000000_255000000:2p:240000000-255000000"
    
#     # "CHR3_115000000_135000000:3p:115000000-135000000"
#     # "CHR3_115000000_135000000:3q:115000000-135000000"
#     # "CHR4_305000000_310000000:4:305000000-310000000"
#     # "CHR5_240000000_250000000:5:240000000-250000000"
#     # "CHR6_135000000_140000000:6:135000000-140000000"
#     # "CHR6_230000000_240000000:6:230000000-240000000"
# )

# # ── Collect CRAM files ────────────────────────────────────────────────────────
# # Accepts either:
# #   - a file containing one CRAM path per line (output of find_cram_files.sh)
# #   - a directory containing CRAM files directly
# if [[ -f "$DATADIR" ]]; then
#     CRAM_FILES=()
#     while IFS= read -r line; do
#         CRAM_FILES+=("$line")
#     done < "$DATADIR"
# elif [[ -d "$DATADIR" ]]; then
#     shopt -s nullglob
#     CRAM_FILES=("$DATADIR"/*.cram)
#     shopt -u nullglob
# else
#     echo "ERROR: '$DATADIR' is neither a file nor a directory."
#     echo "Usage: $0 <cram_paths_file_or_directory> [output_base_dir]"
#     exit 1
# fi

# if [[ ${#CRAM_FILES[@]} -eq 0 ]]; then
#     echo "No CRAM files found via '$DATADIR'"
#     exit 1
# fi

# echo "Found ${#CRAM_FILES[@]} CRAM file(s)."
# echo "Output base directory: $RESDIR"
# echo ""

# # ── Main loop ─────────────────────────────────────────────────────────────────
# for region_entry in "${REGIONS[@]}"; do
#     LABEL="${region_entry%%:*}"          # e.g. CHR1_240000000_255000000
#     REGION="${region_entry#*:}"          # e.g. 1:240000000-255000000

#     OUT_DIR="${RESDIR}/${LABEL}"
#     mkdir -p "$OUT_DIR"

#     echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
#     echo "Region: $REGION  →  $OUT_DIR"
#     echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

#     for cram in "${CRAM_FILES[@]}"; do
#         sample=$(basename "$cram" .cram)
#         out_cram="${OUT_DIR}/${sample}_${LABEL}.cram"


#         echo "  Processing: $sample"

#         -C requires fasta input. Transferred from farm to local dir.
#         samtools view \
#             --reference "$REF" \
#             -C \
#             -o "$out_cram" \
#             "$cram" \
#             "$REGION"
#         out="${OUT_DIR}/${sample}_${LABEL}_coverage.txt"
#         samtools coverage "$cram" -o "$out" --reference "$REF" --region "$REGION"
        
#         out="${OUT_DIR}/${sample}_${LABEL}_depthPlot.txt"
#         samtools coverage "$cram" -o "$out" --reference "$REF" --region "$REGION" --plot-depth

#         if [[ $? -eq 0 ]]; then
#             echo "    Indexing:  $out_cram"
#             samtools index "$out_cram"
#             echo "    Done:      $out_cram"
#         else
#             echo "    ERROR processing $sample for region $REGION"
#         fi
#     done
# done

# echo "All regions complete. Extracts saved under '$RESDIR'."


# ── Regions of interest ───────────────────────────────────────────────────────
# Format: "LABEL:ARM:START-END"
#   LABEL  – output directory / file prefix
#   ARM    – chromosome name as it appears in the CRAM header (e.g. 1p, 1q, 2p)
#   START-END – coordinates in the *original* (unsplit, absolute) space.
#               For q-arm entries the script subtracts the p-arm length
#               automatically, so you always write absolute genomic positions.
REGIONS=(
    # Label is absolute position. Specify whether region falls on p or q arm to adjust coordinates. 
    "CHR1_530000000_600000000:1q:530000000-600000000"
    "CHR1_240000000_255000000:1p:240000000-255000000"
    "CHR2_240000000_255000000:2p:240000000-255000000"
    "CHR3_115000000_135000000:3p:115000000-135000000"
    "CHR4_305000000_310000000:4:305000000-310000000"
    "CHR5_240000000_250000000:5:240000000-250000000"
    "CHR6_135000000_140000000:6:135000000-140000000"
    "CHR6_230000000_240000000:6:230000000-240000000"
)

# ── Collect CRAM files ────────────────────────────────────────────────────────
if [[ -f "$DATADIR" ]]; then
    CRAM_FILES=()
    while IFS= read -r line; do
        [[ -n "$line" ]] && CRAM_FILES+=("$line")
    done < "$DATADIR"
elif [[ -d "$DATADIR" ]]; then
    shopt -s nullglob
    CRAM_FILES=("$DATADIR"/*.cram)
    shopt -u nullglob
else
    echo "ERROR: '$DATADIR' is neither a file nor a directory."
    echo "Usage: $0 <cram_paths_file_or_directory> [output_base_dir]"
    exit 1
fi

if [[ ${#CRAM_FILES[@]} -eq 0 ]]; then
    echo "No CRAM files found via '$DATADIR'"
    exit 1
fi

# ── Build chromosome-length lookup from the first CRAM header ─────────────────
# Saves a TSV  <chrom_name> <TAB> <length>  to $RESDIR/chrom_lengths.tsv
# so it can be inspected / reused between runs.
mkdir -p "$RESDIR"
CHROM_LENGTHS_FILE="${RESDIR}/chrom_lengths.tsv"

echo "Extracting chromosome lengths from header of: ${CRAM_FILES[0]}"
samtools view -H "${CRAM_FILES[0]}" \
    | awk '/^@SQ/ {
        name=""; len=0
        for (i=2; i<=NF; i++) {
            if ($i ~ /^SN:/) name=substr($i,4)
            if ($i ~ /^LN:/) len=substr($i,4)
        }
        if (name != "") print name "\t" len
      }' \
    > "$CHROM_LENGTHS_FILE"

echo "Chromosome lengths saved to: $CHROM_LENGTHS_FILE"
echo ""

# ── Helper: look up a chromosome length from the TSV ─────────────────────────
chrom_length() {
    local chrom="$1"
    awk -v c="$chrom" '$1==c {print $2; exit}' "$CHROM_LENGTHS_FILE"
}

# ── Helper: compute offset for a q-arm coordinate ────────────────────────────
# If the ARM ends in 'q', subtract the matching p-arm length from the
# absolute coordinates so they become p-arm-relative.
# Returns three values via global variables: RESOLVED_ARM, ADJ_START, ADJ_END
resolve_region() {
    local arm="$1"    # e.g. 1q
    local start="$2"  # absolute genomic start
    local end="$3"    # absolute genomic end

    RESOLVED_ARM="$arm"
    ADJ_START="$start"
    ADJ_END="$end"

    if [[ "$arm" =~ ^(.+)q$ ]]; then
        local base_chr="${BASH_REMATCH[1]}"   # e.g. "1"
        local p_arm="${base_chr}p"
        local p_len
        p_len=$(chrom_length "$p_arm")

        if [[ -z "$p_len" ]]; then
            echo "  WARNING: Could not find length for '${p_arm}' in header; using raw coordinates." >&2
            return
        fi

        ADJ_START=$(( start - p_len ))
        ADJ_END=$(( end   - p_len ))

        if (( ADJ_START < 1 )); then
            echo "  WARNING: Adjusted start for ${arm}:${start}-${end} is ${ADJ_START} (< 1)." \
                 "Check that your absolute coordinates are correct." >&2
            ADJ_START=1
        fi

        echo "  Offset: ${arm} absolute ${start}-${end}" \
             "→ ${arm} relative $(( start - p_len ))-$(( end - p_len ))" \
             " (p-arm length = ${p_len})"
    fi
}

# ── Main loop ─────────────────────────────────────────────────────────────────
echo "Found ${#CRAM_FILES[@]} CRAM file(s)."
echo "Output base directory: $RESDIR"
echo ""

for region_entry in "${REGIONS[@]}"; do
    LABEL="${region_entry%%:*}"              # CHR1q_530000000_600000000
    rest="${region_entry#*:}"               # 1q:530000000-600000000
    ARM="${rest%%:*}"                        # 1q
    COORDS="${rest#*:}"                      # 530000000-600000000
    RAW_START="${COORDS%-*}"                 # 530000000
    RAW_END="${COORDS#*-}"                   # 600000000

    # Resolve coordinates (applies p-arm offset for q-arm entries)
    resolve_region "$ARM" "$RAW_START" "$RAW_END"
    # RESOLVED_ARM, ADJ_START, ADJ_END are now set

    REGION_STR="${RESOLVED_ARM}:${ADJ_START}-${ADJ_END}"
    OUT_DIR="${RESDIR}/${LABEL}"
    mkdir -p "$OUT_DIR"

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Region: ${ARM}:${RAW_START}-${RAW_END}  →  query: ${REGION_STR}"
    echo "Output: $OUT_DIR"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    for cram in "${CRAM_FILES[@]}"; do
        sample=$(basename "$cram" .cram)

        out_cram="${OUT_DIR}/${sample}_${LABEL}.cram"

        echo "  Processing: $sample"

        samtools view \
            --reference "$REF" \
            -C \
            -o "$out_cram" \
            "$cram" \
            "$REGION_STR"

        samtools index "$out_cram"

        out="${OUT_DIR}/${sample}_${LABEL}_coverage.txt"
        samtools coverage "$cram" -o "$out" --reference "$REF" --region "$REGION_STR"
        
        out="${OUT_DIR}/${sample}_${LABEL}_depthPlot.txt"
        samtools coverage "$cram" -o "$out" --reference "$REF" --region "$REGION_STR" --plot-depth

        if [[ $? -eq 0 ]]; then
            echo "    Done:      $out_cram"
        else
            echo "    ERROR processing $sample for region $REGION"
        fi

    done
done

echo ""
echo "All regions complete. Results saved under '${RESDIR}'."
