# Scripts for extra experiment and plot

- `false_discovery_causal.R` Analyze the false discoveries in the causal simulation for all non-causal cell types across traits.

- `label_perturb_random_seeds.R` generates random cell id seeds for label perturbation experiments. The sampled cell IDs are used for downstream simulation of label reassignment. The output file also include a data frame with all the details of the randomization (the target cell type names, the number of cells in the target cell types, the output file header, etc).

- `label_perturb_score_calc.R` Calculate the seismic specificity score and other metrics (including Bryios gene specificity, specificity index and DE score) for the label perturbation experiment. The script uses the random cell id seeds generated in `label_perturb_random_seeds.R` to reassign cells as the target cell types. It then calculates the scores for each gene for the target cell types and saves the results in a .tsv file.

- `label_perturb_res.R` Summarise the results from the label perturbation experiment. The script reads the .tsv files generated in `label_perturb_score_calc.R` and calculates the normalized L1 distance between the original scores and the score after label perturbation, for all metrics. The output is Supplementary Figure 1.

- `multi_individuals_analysis.R` Run seismic trait-cell type associationnt analysis for each individual donor in the Tabula muris dataset. The results are in Supplementary Figure 22.

- `pancreatic_score_plot.R` Plot the seismic specificity score for the pancreatic cell types in the TM dataset (Supplementary Figure 2).

- `Saunders_markers_plot.R` Plot the distribution of the seismic specificity score for the subtypes of substantia nigra neurons and the score when merging the three subtypes into one (Supplementary Figure 24). 

- `Saunders_multiple_filter.R` Anlyze the Saunders et al dataset using Parkinson's disease GWAS using three different threshold of cell count filter (10, 20, 50). The results are in Supplementary Figure 23.

- `spearman_analysis.R` Analyze the TM dataset using the 27 traits when specifying a nonparametric Spearman's correlation model. The results are displayed in Supplementary Figure 29.
