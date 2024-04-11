#!/bin/bash

# Check for the right number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 CONTAINER_NAME DATA_PATH OUTPUT_PATH"
    exit 1
fi

CONTAINER_NAME="$1"
DATA_PATH="$2"
OUTPUT_PATH="$3"
IMAGE_NAME="seismic_rt"

# Check if the data path exists
if [ ! -f "$DATA_PATH" ]; then
    echo "Error: Data path '$DATA_PATH' does not exist."
    exit 1
fi

# Check if the output path exists and create if not
if [ ! -d "$OUTPUT_PATH" ]; then
    mkdir -p "$OUTPUT_PATH"
fi

# Initiate container
podman run -it -d --name $CONTAINER_NAME --cpu-quota 100 $IMAGE_NAME

# Copy the data file
data_name=`basename $DATA_PATH`
podman cp "$DATA_PATH" $CONTAINER_NAME:/grain/ql29/podman_file/data/tmp/$data_name

# Run the script
res_name=`basename $OUTPUT_PATH`
podman exec $CONTAINER_NAME Rscript /grain/ql29/podman_file/scripts/ours_running_time.R /grain/ql29/podman_file/data/tmp/$data_name /grain/ql29/podman_file/data/zscore_5 /grain/ql29/podman_file/data/tmp/$res_name

# Copy data to output path
podman cp $CONTAINER_NAME:/grain/ql29/podman_file/data/tmp/$res_name "$OUTPUT_PATH"

# Remove container
podman stop $CONTAINER_NAME
podman rm $CONTAINER_NAME  # Fixed typo in variable name

echo "Script finished successfully."