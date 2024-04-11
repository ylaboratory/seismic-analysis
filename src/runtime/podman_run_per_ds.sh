#!/bin/bash

# Check for the right number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 CONTAINER_NAME DATA_PATH_RDA DATA_PATH_H5AD OUTPUT_PATH"
    exit 1
fi

#get the base path and set it as the working directory
BASE_PATH=$(dirname $(dirname $(realpath $0))) 

CONTAINER_NAME="$1"
DATA_PATH_RDA="$2" #data path of the rdata
DATA_PATH_H5AD="$3" #data path of the h5ad
OUTPUT_PATH="$4" #output path
IMAGE_NAME="seismic_rt"

# Check if the data path exists
if [ ! -f "$DATA_PATH_RDA" ] || [ ! -f "$DATA_PATH_H5AD" ] ; then
    echo "Error: Data '$DATA_PATH_RDA' or "$DATA_PATH_H5AD" do not exist."
    exit 1
fi

# Check if the output path exists and create if not
if [ ! -d "$OUTPUT_PATH" ]; then
    mkdir -p "$OUTPUT_PATH"
fi

# Initiate container
podman run -it -d --name $CONTAINER_NAME --cpu-quota 100 $IMAGE_NAME

#data name: data file name header
data_name_rda=$(basename $DATA_PATH_RDA)
data_name_h5ad=$(basename $DATA_PATH_H5AD)
#data header: data file directory name
data_header_rda=$(basename $(dirname $DATA_PATH_RDA))
data_header_h5ad=$(basename $(dirname $DATA_PATH_H5AD))

#cp data to container
podman cp "$BASE_PATH"/"$DATA_PATH_RDA" $CONTAINER_NAME:data/tmp/"$data_header_rda"."$data_name_rda"
podman cp "$BASE_PATH"/"$DATA_PATH_H5AD" $CONTAINER_NAME:data/tmp/"$data_header_h5ad"."$data_name_h5ad"

# Run the script
res_name=$(basename $OUTPUT_PATH)
podman exec "$CONTAINER_NAME" Rscript scripts/magma_data_gen.R data/tmp/"$data_header_rda"."$data_name_rda" data/raw_5 data/tmp/"$res_name".magma.rda data/tmp/"$res_name"_test 
podman exec "$CONTAINER_NAME" Rscript scripts/fuma_data_gen.R data/tmp/"$data_header_rda"."$data_name_rda" data/raw_5 data/tmp/"$res_name".fuma.rda data/tmp/"$res_name"_test 
podman exec "$CONTAINER_NAME" Rscript scripts/seismic_runtime.R data/tmp/"$data_header_rda"."$data_name_rda" data/zscore_5 data/tmp/"$res_name".seismic.rda 
podman exec "$CONTAINER_NAME" python scripts/scdrs_running_time.py data/tmp/"$data_header_h5ad"."$data_name_h5ad" data/scdrs_5 data/tmp/"$res_name".scdrs.joblib

# Copy data to output path
podman cp $CONTAINER_NAME:/grain/ql29/podman_file/data/tmp/$res_name".magma.rda" "$BASE_PATH"/"$OUTPUT_PATH"
podman cp $CONTAINER_NAME:/grain/ql29/podman_file/data/tmp/$res_name".fuma.rda" "$BASE_PATH"/"$OUTPUT_PATH"
podman cp $CONTAINER_NAME:/grain/ql29/podman_file/data/tmp/$res_name".seismic.rda" "$BASE_PATH"/"$OUTPUT_PATH"
podman cp $CONTAINER_NAME:/grain/ql29/podman_file/data/tmp/$res_name".scdrs.joblib" "$BASE_PATH"/"$OUTPUT_PATH"

# Remove container
podman stop $CONTAINER_NAME
podman rm $CONTAINER_NAME  # Fixed typo in variable name