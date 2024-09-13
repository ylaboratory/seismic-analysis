#!/bin/bash

#set number of cores
NUM_CORE=15

#standard 
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/standard/perturbed_expr/parameter_df.txt results/causal_sim/standard es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/standard/perturbed_expr/parameter_df.txt \
results/causal_sim/standard es_name $NUM_CORE > results/causal_sim/standard/scdrs_process.log &

#unround
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/standard/perturbed_expr_unround/parameter_df.txt results/causal_sim/standard_unround es_name TRUE $NUM_CORE &


#vary gene ratio
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/varying_gene_ratio/perturbed_expr/parameter_df.txt results/causal_sim/varying_gene_ratio gr_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/varying_gene_ratio/perturbed_expr/parameter_df.txt \
results/causal_sim/varying_gene_ratio gr_name $NUM_CORE > results/causal_sim/varying_gene_ratio/scdrs_process.log &

#vary cell ratio
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/varying_cell_ratio/perturbed_expr/parameter_df.txt results/causal_sim/varying_cell_ratio cr_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/varying_cell_ratio/perturbed_expr/parameter_df.txt \
results/causal_sim/varying_cell_ratio cr_name $NUM_CORE > results/causal_sim/varying_cell_ratio/scdrs_process.log &

#vary cell ratio - only for disease genes
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/varying_cell_ratio/perturbed_expr_all_cells/parameter_df.txt results/causal_sim/varying_cell_ratio_all_cells cr_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/varying_cell_ratio/perturbed_expr_all_cells/parameter_df.txt \
results/causal_sim/varying_cell_ratio_all_cells cr_name $NUM_CORE > results/causal_sim/varying_cell_ratio_all_cells/scdrs_process.log &

#multiple cell types

#correlated cell types

#multiple associated cell types

##B cells (real cell types)


###other scenario
#scDesign3
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_standard/perturbed_expr/parameter_df.txt results/causal_sim/sc3_standard es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_standard/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_standard es_name $NUM_CORE > results/causal_sim/sc3_standard/sc3_scdrs_process.log &

#scDesign3 B cells
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_ct/perturbed_expr/parameter_df.txt results/causal_sim/sc3_ct es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_ct/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_ct es_name $NUM_CORE > results/causal_sim/sc3_ct/sc3_scdrs_process.log &

#scDesign3 neuron
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_neuro/perturbed_expr/parameter_df.txt results/causal_sim/sc3_neuro es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_neuro/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_neuro es_name $NUM_CORE > results/causal_sim/sc3_neuro/sc3_scdrs_process.log &

#scDesign3: different numbers of genes
##100 causal genes
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_gene_num/gene_100/perturbed_expr/parameter_df.txt results/causal_sim/sc3_gene_num/gene_100 es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_gene_num/gene_100/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_gene_num/gene_100 es_name $NUM_CORE > results/causal_sim/sc3_gene_num/gene_100/sc3_scdrs_process.log &

##400 causal genes
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_gene_num/gene_400/perturbed_expr/parameter_df.txt results/causal_sim/sc3_gene_num/gene_400 es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_gene_num/gene_400/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_gene_num/gene_400 es_name $NUM_CORE > results/causal_sim/sc3_gene_num/gene_400/sc3_scdrs_process.log &

#scDesign3: different ratio of cells
##20 cells
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_cell_num/cell_20/perturbed_expr/parameter_df.txt results/causal_sim/sc3_cell_num/cell_20 es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_cell_num/cell_20/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_cell_num/cell_20 es_name $NUM_CORE > results/causal_sim/sc3_cell_num/cell_20/sc3_scdrs_process.log &

##50 cells
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_cell_num/cell_50/perturbed_expr/parameter_df.txt results/causal_sim/sc3_cell_num/cell_50 es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_cell_num/cell_50/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_cell_num/cell_50 es_name $NUM_CORE > results/causal_sim/sc3_cell_num/cell_50/sc3_scdrs_process.log &

# scDesign3: probabilistic gene sets
## 200 genes
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_prob/gene_200/perturbed_expr/parameter_df.txt results/causal_sim/sc3_prob/gene_200 es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_prob/gene_200/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_prob/gene_200 es_name $NUM_CORE > results/causal_sim/sc3_prob/gene_200/sc3_scdrs_process.log &

## 300 genes
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_prob/gene_300/perturbed_expr/parameter_df.txt results/causal_sim/sc3_prob/gene_300 es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_prob/gene_300/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_prob/gene_300 es_name $NUM_CORE > results/causal_sim/sc3_prob/gene_300/sc3_scdrs_process.log &

## 400 genes
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_prob/gene_400/perturbed_expr/parameter_df.txt results/causal_sim/sc3_prob/gene_400 es_name TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_prob/gene_400/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_prob/gene_400 es_name $NUM_CORE > results/causal_sim/sc3_prob/gene_400/sc3_scdrs_process.log &


#scDesign3 expanded
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_prob_extended/tot_1000/perturbed_expr/parameter_df.txt results/causal_sim/sc3_prob_extended/tot_1000 output_header FALSE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_prob_extended/tot_1000/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_prob_extended/tot_1000 output_header FALSE $NUM_CORE > results/causal_sim/sc3_prob_extended/tot_1000/sc3_scdrs_process.log &

#2000
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_prob_extended/tot_2000/perturbed_expr/parameter_df.txt results/causal_sim/sc3_prob_extended/tot_2000 output_header FALSE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_prob_extended/tot_2000/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_prob_extended/tot_2000 output_header FALSE $NUM_CORE > results/causal_sim/sc3_prob_extended/tot_2000/sc3_scdrs_process.log &


#500
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_prob_extended/tot_500/perturbed_expr/parameter_df.txt results/causal_sim/sc3_prob_extended/tot_500 output_header TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_prob_extended/tot_500/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_prob_extended/tot_500 output_header FALSE $NUM_CORE > results/causal_sim/sc3_prob_extended/tot_500/sc3_scdrs_process.log &


#5000
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_prob_extended/tot_5000/perturbed_expr/parameter_df.txt results/causal_sim/sc3_prob_extended/tot_5000 output_header TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_prob_extended/tot_5000/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_prob_extended/tot_5000 output_header FALSE target_cell_type $NUM_CORE > results/causal_sim/sc3_prob_extended/tot_5000/sc3_scdrs_process.log &


#new
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_gene_num_new/gene_400/perturbed_expr/parameter_df.txt results/causal_sim/sc3_gene_num_new/gene_400 es_name TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_gene_num_new/gene_400/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_gene_num_new/gene_400 es_name $NUM_CORE > results/causal_sim/sc3_gene_num_new/gene_400/sc3_scdrs_process.log &


##multiple cell types
###500 genes
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_multi_ct_500/random/perturbed_expr/parameter_df.txt results/causal_sim/sc3_multi_ct_500/random output_header target_cell_type TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_multi_ct_500/random/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/random output_header FALSE $NUM_CORE > results/causal_sim/sc3_multi_ct_500/random/sc3_scdrs_process.log &

#no overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_multi_ct_500/no_overlap/perturbed_expr/parameter_df.txt results/causal_sim/sc3_multi_ct_500/no_overlap output_header target_cell_type TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_multi_ct_500/no_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/no_overlap output_header FALSE $NUM_CORE > results/causal_sim/sc3_multi_ct_500/no_overlap/sc3_scdrs_process.log &

#ct overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_multi_ct_500/ct_overlap/perturbed_expr/parameter_df.txt results/causal_sim/sc3_multi_ct_500/ct_overlap output_header target_cell_type TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_multi_ct_500/ct_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/ct_overlap output_header FALSE $NUM_CORE > results/causal_sim/sc3_multi_ct_500/ct_overlap/sc3_scdrs_process.log &

#causal overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_multi_ct_500/causal_overlap/perturbed_expr/parameter_df.txt results/causal_sim/sc3_multi_ct_500/causal_overlap output_header target_cell_type TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_multi_ct_500/causal_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/causal_overlap output_header FALSE $NUM_CORE > results/causal_sim/sc3_multi_ct_500/causal_overlap/sc3_scdrs_process.log &


###10 cell types
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_10_ct/random/perturbed_expr/parameter_df.txt results/causal_sim/sc3_10_ct/random output_header target_cell_type TRUE FALSE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_10_ct/random/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_10_ct/random output_header target_cell_type TRUE $NUM_CORE > results/causal_sim/sc3_10_ct/random/sc3_scdrs_process.log &

#no overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_10_ct/no_overlap/perturbed_expr/parameter_df.txt results/causal_sim/sc3_10_ct/no_overlap output_header target_cell_type TRUE FALSE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_10_ct/no_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_10_ct/no_overlap output_header target_cell_type TRUE $NUM_CORE > results/causal_sim/sc3_10_ct/no_overlap/sc3_scdrs_process.log &

#ct overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_10_ct/ct_overlap/perturbed_expr/parameter_df.txt results/causal_sim/sc3_10_ct/ct_overlap output_header target_cell_type TRUE FALSE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_10_ct/ct_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_10_ct/ct_overlap output_header target_cell_type TRUE $NUM_CORE > results/causal_sim/sc3_10_ct/ct_overlap/sc3_scdrs_process.log &

#causal overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_10_ct/causal_overlap/perturbed_expr/parameter_df.txt results/causal_sim/sc3_10_ct/causal_overlap output_header target_cell_type TRUE FALSE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_10_ct/causal_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_10_ct/causal_overlap output_header target_cell_type TRUE $NUM_CORE > results/causal_sim/sc3_10_ct/causal_overlap/sc3_scdrs_process.log &

#large effect

#limited number of genes with large effect


##get the scDesign3 model
for i in {2..10}
do
nohup Rscript src/data-cleaning/fit_scDesign3.R data/expr/null_sim/expr_rda_rs/expr_ds_${i}.rda data/expr/causal_sim/standard_scdesign3_model/expr_ds_${i}.cell_anno.txt \
  data/expr/causal_sim/standard_scdesign3_model/expr_ds_${i}.scdesign3_para.rda 10 &
done
#

for i in {2..10}
do
nohup Rscript src/data-cleaning/fit_scDesign3.R data/expr/null_sim/expr_rda_rs/expr_ds_${i}.rda \
  data/expr/causal_sim/ct_scdesign3_model/expr_ds_${i}.scdesign3_para.rda 10 &
done