#!/bin/bash
##get the scDesign3 model
for i in {1..10}
do
nohup Rscript src/data-cleaning/fit_scDesign3.R data/expr/null_sim/expr_rda_rs/expr_ds_${i}.rda data/expr/causal_sim/standard_scdesign3_model/expr_ds_${i}.cell_anno.txt \
  data/expr/causal_sim/standard_scdesign3_model/expr_ds_${i}.scdesign3_para.rda 10 &
done


NUM_CORE=15
#scDesign3 expanded
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_prob_extended/tot_1000/perturbed_expr/parameter_df.txt\
results/causal_sim/sc3_prob_extended/tot_1000 output_header target_cell_type FALSE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_prob_extended/tot_1000/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_prob_extended/tot_1000 output_header target_cell_type FALSE $NUM_CORE > results/causal_sim/sc3_prob_extended/tot_1000/sc3_scdrs_process.log &

##multiple cell types
#random 
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_multi_ct_500/random/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/random output_header target_cell_type TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_multi_ct_500/random/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/random output_header target_cell_type FALSE $NUM_CORE > results/causal_sim/sc3_multi_ct_500/random/sc3_scdrs_process.log &

#no overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_multi_ct_500/no_overlap/perturbed_expr/parameter_df.txt\
results/causal_sim/sc3_multi_ct_500/no_overlap output_header target_cell_type TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_multi_ct_500/no_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/no_overlap output_header target_cell_type FALSE $NUM_CORE > results/causal_sim/sc3_multi_ct_500/no_overlap/sc3_scdrs_process.log &

#ct overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_multi_ct_500/ct_overlap/perturbed_expr/parameter_df.txt\
results/causal_sim/sc3_multi_ct_500/ct_overlap output_header target_cell_type TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_multi_ct_500/ct_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/ct_overlap output_header target_cell_type FALSE $NUM_CORE > results/causal_sim/sc3_multi_ct_500/ct_overlap/sc3_scdrs_process.log &

#causal overlap
nohup Rscript src/causal-sim/causal_sim.R data/expr/causal_sim/sc3_multi_ct_500/causal_overlap/perturbed_expr/parameter_df.txt\
results/causal_sim/sc3_multi_ct_500/causal_overlap output_header target_cell_type TRUE TRUE $NUM_CORE &

nohup python -u src/causal-sim/scdrs_causal_sim.py data/expr/causal_sim/sc3_multi_ct_500/causal_overlap/perturbed_expr/parameter_df.txt \
results/causal_sim/sc3_multi_ct_500/causal_overlap output_header target_cell_type FALSE $NUM_CORE > results/causal_sim/sc3_multi_ct_500/causal_overlap/sc3_scdrs_process.log &

