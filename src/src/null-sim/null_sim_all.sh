#!/bin/bash

#set number of cores
NUM_CORE=20

#run for scrambled expression data
#run cross every sample and every seed: fuma magma and seismic
for i in {1..10}
do
    nohup Rscript src/null-sim/seismic_fuma_magma_null_sim.R data/expr/null_sim/expr_rda_rs/expr_ds_"$i".rda \
        data/expr/null_sim/seed_table/sample_cell_idx."$i".txt data/gwas/null_sim/zscore/gs_"$i".genes.out \
        data/gwas/null_sim/magma_raw/gs_"$i".genes.raw data/log/null_sim/ds_"$i".all.process.log  "$NUM_CORE"\
        results/null_sim/expr_rs/new_ds_"$i".null_res.txt data/temp data/ref/mapping/mmu_hsa_mapping.rda &
done

for i in {1..10}
do
    nohup Rscript src/null-sim/fuma_null_sim.R data/expr/null_sim/expr_rda_rs/expr_ds_"$i".rda \
        data/expr/null_sim/seed_table/sample_cell_idx."$i".txt data/gwas/null_sim/zscore/gs_"$i".genes.out \
        data/gwas/null_sim/magma_raw/gs_"$i".genes.raw data/log/null_sim/ds_"$i".fuma.process.log  "$NUM_CORE"\
        results/null_sim/expr_rs/fuma_ds_"$i".null_res.txt data/temp  &
done


for i in {1..10}
do
    nohup Rscript src/null-sim/magma_null_sim.R data/expr/null_sim/expr_rda_rs/expr_ds_"$i".rda \
        data/expr/null_sim/seed_table/sample_cell_idx."$i".txt data/gwas/null_sim/zscore/gs_"$i".genes.out \
        data/gwas/null_sim/magma_raw/gs_"$i".genes.raw data/log/null_sim/ds_"$i".magma.process.log  "$NUM_CORE"\
        results/null_sim/expr_rs/magma_ds_"$i".null_res.txt data/temp data/ref/mapping/mmu_hsa_mapping.rda &
done

#run cross every sample and every seed: scdrs
for i in {3..10}
do
    nohup python -u src/null-sim/scdrs_null_sim.py data/expr/null_sim/expr_h5ad_rs/expr_ds_"$i".h5ad \
        data/expr/null_sim/seed_table/sample_cell_idx."$i".txt data/gwas/null_sim/scdrs_gs/gs_"$i".gs "$NUM_CORE" \
        results/null_sim/expr_rs/new_ds_"$i".scdrs_null_res.txt > data/log/null_sim/ds_"$i".scdrs.process.log &
done


#run for scrambled gene set
for i in {1..10}
do
   nohup Rscript src/null-sim/seismic_null_sim.R data/expr/null_sim/expr_rda/expr_ds_"$i".rda \
        data/expr/null_sim/seed_table/sample_cell_idx."$i".txt  data/gwas/null_sim/zscore_rs/gs_"$i".genes.out \
        data/log/null_sim/gsrs_ds_"$i".all.process.log   "$NUM_CORE" \
        results/null_sim/gs_rs/new_ds_"$i".null_res.txt &
done

#run cross every sample and every seed: scdrs
for i in {1..10}
do
    nohup python -u src/null_sim/scdrs_null_sim.py data/expr/null_sim/expr_h5ad/expr_ds_"$i".h5ad \
        data/expr/null_sim/seed_table/sample_cell_idx."$i".txt data/gwas/null_sim/scdrs_gs_rs/gs_"$i".gs "$NUM_CORE" \
        results/null_sim/gs_rs/new_ds_"$i".scdrs_null_res.txt > data/log/null_sim/gsrs_ds_"$i".scdrs.process.log
done