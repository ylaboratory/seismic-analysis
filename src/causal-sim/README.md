The directory `src/causal-sim` is used for simulating cell groups trait-specific gene expression and assessing the performance of the different methods in detecting them. All synthetic data are generated using the [scDesign3](https://www.nature.com/articles/s41587-023-01772-1) approach. The step-by-step usage for generating the synthetic data works as follows:

- Sample cells from datasets as the cell group to be perturbed. This is achieved in the `sample_cells.R` script. After running the script, the new cell type information will be printed and saved. 

- Fit corresponding scDeisng3 models for datasets with new cell type information. This is achieved in the `fit_scDesign3.R` script. The command line arguments for the script are the dataset file path and the cell type information file path and the final output directory.  After fitting the 10 scDesign3 models, each individual model will be merged into a whole object with `merge_scDesign3_model_obj.R`.

- Simulate new expression profiles with the fitted scDesign3 models. In the `causal_sim_data_gene.R` script, the target gene sets (genes to be perturbed) are sampled based on different sampling strategies. After determing the gene sets and the target cell lists, the saved scDesign3 models are reloaded and set as the parameters for the simulation. The synthetic data matices and the target gene set information are saved.

- Run the methods on the simulated data and asses the performance. `causal_sim.R` are used to run across  *seismic*, S-MAGMA, and FUMA in batch for all the simulated datasets, while `scdsrs_causal_sim.py` is used to run the scDSRS method. The two scripts are scheduled as parallel program to accelerate the process. The output files include reported results for all cell types plus *seismic*'s influential gene analysis for the target cell type.

- The shell script `causal_sim.sh` contains all commands for fitting the scDesign3 models (run `fit_scDesign3`, `causal_sim.R` and `scdrs_causal_sim.py`) and assessing the performance. 

