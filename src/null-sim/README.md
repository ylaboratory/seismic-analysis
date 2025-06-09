# Null Simulation for Framework Calibration Assessment

This directory contains scripts for performing null simulations to assess the calibration of different cell type-trait association frameworks.
The simulation process uses randomly generated seeds, shuffled expression data, and random trait MAGMA gene-level z-scores to calculate P-values for the null hypothesis.
The process is repeated 10,000 times to generate a null distribution of P-values. The scripts are organized as follows:

- `null_sim_all_frameworks.sh`: A shell script containing commands to generate null simulation results for the scDRS, FUMA, and S-MAGMA frameworks.

- `scdrs_null_sim.py`: A Python script for running the null simulation for the scDRS framework. It:
  - Reads shuffled expression data and random trait scDRS gene set files
  - Calculates P-values for the null hypothesis
  - Saves results in the `results/null_sim` directory

- `seismic_fuma_magma_null_sim.R`: An R script for running null simulations for the FUMA and S-MAGMA frameworks, as well as the seismic analysis. It:
  - Reads shuffled expression data
  - Uses random seed samples as target cell types
  - Generates corresponding S-MAGMA and FUMA analysis files
  - Performs _seismic_ analysis
  - Saves results in the `results/null_sim` directory

For each simulation run:
1. A random seed is extracted, representing a sample of 100 cells as the target cell type
2. Shuffled expression data and random trait MAGMA gene-level z-scores are used
3. P-values are calculated for the null hypothesis
4. Results are saved for further analysis and comparison