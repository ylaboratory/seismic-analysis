#!/bin/bash

#run cross every sample size and every data set
#sample size
sample_size=(10k 25k 50k 100k 150k 200k 250k 300k)
data_set=(ds_1 ds_2 ds_3 ds_4 ds_5)

for i in "${sample_size[@]}";
do
    for j in "${data_set[@]}";
    do
        echo "Running for sample size $i and data set $j"
        #if the output path does not exist, create it
        if [ ! -d results/runtime/"$i"/"$j" ]; then
            mkdir -p results/runtime/"$i"/"$j"
        fi
        #run the script
        nohup bash src/runtime/podman_run_per_ds.sh container_"$i"_"$j"  data/expr/runtime/expr_rda/"$i"/sample."$j".rda data/expr/runtime/expr_h5ad/"$i"/sample."$j".h5ad results/runtime/"$i"/"$j" > data/log/runtime/"$i"_"$j".log 2>&1 &
        echo "Finished running for sample size $i and data set $j"
    done
done    

