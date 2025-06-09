#!/bin/bash

print_help() {
    echo "
Usage: $0 -l SNP_LOC_FILE -p SNP_P_VALUE_FILE -o OUTPUT_DIR -s COHORT_SIZE -m MAGMA_PATH -g GENE_LOC_FILE -b B_FILE [-w WINDOW_SIZES] [-h]

This script performs MAGMA annotation and analysis on SNP data.

Required parameters:
    -l | SNP location file
    -p | SNP p-value file
    -o | Output directory
    -s | Cohort size (use 'cX' format to specify column X containing sample size)
    -m | Path to MAGMA binary
    -g | Path to gene location file
    -b | Path to g1000_eur files (without extension)

Optional parameters:
    -w | Window sizes (default: 35,10). Use colon-separated list for multiple sizes, e.g., '35,10:15,20:40,10'
    -h | Display this help message
"
}

# Default values
#MAGMA_PATH="bin/magma/magma"
#GENE_LOC_FILE="data/ref/magma/NCBI37.3.gene.loc"
#B_FILE="data/ref/magma/g1000_eur/g1000_eur"
#WINDOW_SIZES="35,10"


# Default value for window sizes
WINDOW_SIZES="35,10"

# Parse command-line arguments
while getopts "l:p:o:s:m:g:b:w:h" opt; do
    case $opt in
        l) SNP_LOC_FILE="$OPTARG" ;;
        p) SNP_P_VALUE_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) COHORT_SIZE="$OPTARG" ;;
        m) MAGMA_PATH="$OPTARG" ;;
        g) GENE_LOC_FILE="$OPTARG" ;;
        b) B_FILE="$OPTARG" ;;
        w) WINDOW_SIZES="$OPTARG" ;;
        h) print_help; exit 0 ;;
        *) echo "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
    esac
done

# Check for required arguments
if [ -z "$SNP_LOC_FILE" ] || [ -z "$SNP_P_VALUE_FILE" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$COHORT_SIZE" ] || [ -z "$MAGMA_PATH" ] || [ -z "$GENE_LOC_FILE" ] || [ -z "$B_FILE" ]; then
    echo "Error: Missing required arguments" >&2
    print_help
    exit 1
fi

HEADER=$(basename "$SNP_P_VALUE_FILE" | cut -d'.' -f1)
echo "this is header: $HEADER"
echo "${OUTPUT_DIR}/${HEADER}.${ws_fn}"

# Convert window sizes string to array
ws_set=(${WINDOW_SIZES//,/ })

# MAGMA annotation
echo "Starting MAGMA annotation..."
pids=()
for ws in "${ws_set[@]}"; do
    ws_fn=${ws//,/.}
    echo "this is ws: $ws"
    echo "this is fn: $ws_fn"
    echo "this is output file: ${OUTPUT_DIR}/${HEADER}.${ws_fn}"
    echo "this is output directory: $OUTPUT_DIR"
    echo "this is output header: $HEADER"
    "${MAGMA_PATH}" --annotate window="${ws}" --snp-loc "${SNP_LOC_FILE}" --gene-loc "${GENE_LOC_FILE}" --out "${OUTPUT_DIR}/${HEADER}.${ws_fn}" &
    pids+=($!)
done

for pid in "${pids[@]}"; do
    wait "$pid"
done

echo "Finished MAGMA annotation"

# MAGMA analysis
echo "Starting MAGMA analysis..."
analysis_pids=()
for ws in "${ws_set[@]}"; do
    ws_fn=${ws//,/.}
    
    if [[ "$COHORT_SIZE" =~ ^([Cc])([0-9]+)$ ]]; then
        ncol="${BASH_REMATCH[2]}"
        echo "Using column $ncol for sample size"
        total_columns=$(awk '{print NF; exit}' "$SNP_P_VALUE_FILE")
        if [ "$ncol" -gt "$total_columns" ]; then
            echo "Error: COHORT_SIZE value exceeds the total number of columns in SNP_P_VALUE_FILE." >&2
            exit 1
        fi
        "${MAGMA_PATH}" --bfile "${B_FILE}" --gene-annot "${OUTPUT_DIR}/${HEADER}.${ws_fn}.genes.annot" --out "${OUTPUT_DIR}/${HEADER}.${ws_fn}" --pval "${SNP_P_VALUE_FILE}" ncol="$ncol" &
    else
        "${MAGMA_PATH}" --bfile "${B_FILE}" --gene-annot "${OUTPUT_DIR}/${HEADER}.${ws_fn}.genes.annot" --out "${OUTPUT_DIR}/${HEADER}.${ws_fn}" --pval "${SNP_P_VALUE_FILE}" N="${COHORT_SIZE}" &
    fi
    analysis_pids+=($!)
done

for pid in "${analysis_pids[@]}"; do
    wait "$pid"
done

echo "Finished MAGMA analysis"