#!/bin/bash
granularity=("fine_cluster" "subclass" "region_class" "region_subclass" "region_cluster")

##### run FUMA and S-MAGMA for the TM FACS datasets and TM droplet datasets
bash src/tools/FUMA_MAGMA_batch.sh -m MAGMA -e data/expr/Tabula_muris/tm_facs.top10_magma.txt -g data/gwas/tm_gwas/magma_raw -o results/Tabula_muris/FACS/S-MAGMA
bash src/tools/FUMA_MAGMA_batch.sh -m FUMA -e data/expr/Tabula_muris/tm_facs.fuma.txt -g data/gwas/tm_gwas/magma_raw -o results/Tabula_muris/FACS/FUMA
bash src/tools/FUMA_MAGMA_batch.sh -m MAGMA -e data/expr/Tabula_muris/tm_droplet.top10_magma.txt -g data/gwas/tm_gwas/magma_raw -o results/Tabula_muris/droplet/S-MAGMA
bash src/tools/FUMA_MAGMA_batch.sh -m FUMA -e data/expr/Tabula_muris/tm_droplet.fuma.txt -g data/gwas/tm_gwas/magma_raw -o results/Tabula_muris/droplet/FUMA

#run scDRS for the TM data FACS set
conda activate qlrisotto #activate your own environment
all_gs_files=$(ls data/gwas/tm_gwas/scdrs_gs/*.gs)
for f in $all_gs_files; do
    g_prefix=$(basename "$f" | cut -d'.' -f1)
    python bin/scdrs compute-score data/expr/Tabula_muris/facs.clean.h5ad mouse "$f" human results/Tabula_muris/FACS/scDRS --flag-raw-count True
    python bin/scdrs perform-downstream --h5ad-file data/expr/Tabula_muris/facs.clean.h5ad --score-file results/Tabula_muris/FACS/scDRS/"$g_prefix".full_score.gz\
    --group-analysis "cluster_name" --out-folder results/Tabula_muris/FACS/scDRS --flag_raw_count True
done

#for droplet
for f in $all_gs_files; do
    g_prefix=$(basename "$f" | cut -d'.' -f1)
    python bin/scdrs compute-score data/expr/Tabula_muris/droplet.clean.h5ad mouse "$f" human results/Tabula_muris/droplet/scDRS --flag-raw-count True
    python bin/scdrs perform-downstream --h5ad-file data/expr/Tabula_muris/droplet.clean.h5ad --score-file results/Tabula_muris/droplet/scDRS/"$g_prefix".full_score.gz\
    --group-analysis "cluster_name" --out-folder results/Tabula_muris/droplet/scDRS --flag_raw_count True
done

#### run FUMA and S-MAGMA for the TS data sets
bash src/tools/FUMA_MAGMA_batch.sh -m MAGMA -e data/expr/Tabula_sapiens/ts.top10_magma.txt -g data/gwas/tm_gwas/magma_raw -o results/Tabula_sapiens/FACS/S-MAGMA
bash src/tools/FUMA_MAGMA_batch.sh -m FUMA -e data/expr/Tabula_sapiens/ts.fuma.txt -g data/gwas/tm_gwas/magma_raw -o results/Tabula_sapiens/FACS/FUMA

#run scDRS for the TM data FACS set
conda activate qlrisotto #activate your own environment
all_gs_files=$(ls data/gwas/tm_gwas/scdrs_gs/*.gs)
for f in $all_gs_files; do
    g_prefix=$(basename "$f" | cut -d'.' -f1)
    python bin/scdrs compute-score data/expr/Tabula_muris/facs.clean.h5ad mouse "$f" human results/Tabula_muris/FACS/scDRS --flag-raw-count True
    python bin/scdrs perform-downstream --h5ad-file data/expr/Tabula_muris/facs.clean.h5ad --score-file results/Tabula_muris/FACS/scDRS/"$g_prefix".full_score.gz\
    --group-analysis "cluster_name" --out-folder results/Tabula_muris/FACS/scDRS --flag_raw_count True
done


#### Run FUMA and S-MAGMA for the Saunders et al. data set
for g in "${granularity[@]}"; do
    bash src/tools/FUMA_MAGMA_batch.sh -m MAGMA -e data/expr/Saunders/Saunders."$g".top10_magma.txt -g data/gwas/brain_gwas/magma_raw -o results/Saunders/"$g"/S-MAGMA
    bash src/tools/FUMA_MAGMA_batch.sh -m FUMA -e data/expr/Saunders/Saunders."$g".fuma.txt -g data/gwas/brain_gwas/magma_raw -o results/Saunders/"$g"/FUMA
done

for g in "${granularity[@]}"; do
    nohup bash src/tools/FUMA_MAGMA_batch.sh -m MAGMA -e data/expr/Saunders/new_Saunders."$g".top10_magma.txt -g data/gwas/brain_gwas/magma_raw -o results/new_Saunders/"$g"/S-MAGMA &
    nohup bash src/tools/FUMA_MAGMA_batch.sh -m FUMA -e data/expr/Saunders/new_Saunders."$g".fuma.txt -g data/gwas/brain_gwas/magma_raw -o results/new_Saunders/"$g"/FUMA &
done

# Run scDRS for the Saunders et al. data set
conda activate qlrisotto #activate your own environment
#list all files in data/gwas/brain_gwas/scdrs_gs
all_gs_files=$(ls data/gwas/brain_gwas/scdrs_gs/*.gs)
for g in "${granularity[@]}"; do
    #extract the prefix of the gs file
    g_prefix=$(basename "$f" | cut -d'.' -f1)
    for f in $all_gs_files; do
        python bin/scdrs compute-score data/expr/Saunders/Saunders."$g".clean.h5ad mouse "$f" human results/Saunders/"$g"/scDRS --flag-raw-count True
        python bin/scdrs perform-downstream --h5ad-file data/expr/Saunders/Saunders."$g".clean.h5ad --score-file results/Sanders/"$g"/scDRS/"$g_prefix".full_score.gz\
        --group-analysis "$g" --out-folder results/Saunders/"$g"/scDRS --flag_raw_count True
    done
done
conda deactivate 
