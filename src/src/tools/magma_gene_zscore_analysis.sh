#!/bin/bash

# Check for the right number of arguments
if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 SNP_LOC_FILE SNP_P_VALUE_FILE OUTPUT_DIR COHORT_SIZE [SPECIFIC_WS]"
    exit 1
fi

SNP_LOC_FILE="$1"
SNP_P_VALUE_FILE="$2"
OUTPUT_DIR="$3"
COHORT_SIZE="$4"
SPECIFIC_WS="${5:-}"

MAGMA_PATH="bin/magma/magma"
GENE_LOC_FILE="data/ref/magma/NCBI37.3.gene.loc"
B_FILE="data/ref/magma/g1000_eur/g1000_eur"

HEADER=$(basename "$SNP_P_VALUE_FILE" | cut -d'.' -f1)

# Define window size set
ws_set=(
    "100,50"
    "100,10"
    "50,50"
    "50,10"
    "35,10"
    "20,5"
    "10,10"
    "10,1.5"
)

# If a specific window size is given, filter ws_set to include only that window size
if [ -n "$SPECIFIC_WS" ]; then
    ws_set=("$SPECIFIC_WS")
fi

pids=()
for ws in "${ws_set[@]}"; do
    ws_fn=$(echo "$ws" | sed 's/,/./g')
    "${MAGMA_PATH}" --annotate window="${ws}" --snp-loc "${SNP_LOC_FILE}" --gene-loc "${GENE_LOC_FILE}" --out "${OUTPUT_DIR}/${HEADER}.${ws_fn}" &
    pids+=($!)
done

for pid in "${pids[@]}"; do
    wait "$pid"
done

echo "Finished doing MAGMA annotation"

analysis_pids=()
for ws in "${ws_set[@]}"; do
    ws_fn=$(echo "$ws" | sed 's/,/./g')
    
    # if "COHORT_SIZE" starts with C or c, make the later to become the column containing n
    if [[ "$COHORT_SIZE" =~ ^([Cc])([0-9]+)$ ]]; then
        ncol="${BASH_REMATCH[2]}"
        echo "ncol is $ncol"
        # ncol <= total number of columns
        total_columns=$(awk '{print NF}' "$SNP_P_VALUE_FILE" | head -n 1)
        if [ "$ncol" -gt "$total_columns" ]; then
            echo "Error: COHORT_SIZE value exceeds the total number of columns in SNP_P_VALUE_FILE."
            exit 1
        fi
        "${MAGMA_PATH}" --bfile "${B_FILE}" --gene-annot "${OUTPUT_DIR}/${HEADER}.${ws_fn}.genes.annot" --out "${OUTPUT_DIR}/${HEADER}.${ws_fn}" --pval "${SNP_P_VALUE_FILE}" ncol="$ncol" &
        analysis_pids+=($!)
    else
        "${MAGMA_PATH}" --bfile "${B_FILE}" --gene-annot "${OUTPUT_DIR}/${HEADER}.${ws_fn}.genes.annot" --out "${OUTPUT_DIR}/${HEADER}.${ws_fn}" --pval "${SNP_P_VALUE_FILE}" N="${COHORT_SIZE}" &
        analysis_pids+=($!)
    fi
    
done

for pid in "${analysis_pids[@]}"; do
    wait "$pid"
done

echo "Finished doing MAGMA analysis"
