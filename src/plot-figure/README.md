# Generate all figures in the paper using the result files
The `plot-figure` directory contains the scripts for generating all the figures in the paper. The scripts are organized as follows:

- `null_sim_plot.R` processes the null simulation results for the scDRS, FUMA, and S-MAGMA frameworks and generates the calibration plots (Figure 2A, Supplementary Figure 1) for the three frameworks.

- `causal_sim_plot.R` processes the causal simulation results for the four frameworks in various conditions and generates the power plots (Figure 2B, Supplementary Figure 2-3) and the influential gene plots (Figure 2C, Supplementary figure 4) for the four frameworks.

- `runtime_plot.R` processes the runtime evaluation results for the four frameworks and generates the runtime plots (Figure 2D) for the four frameworks. 

- `analysis_TM_FACS_plot.R` processes the analysis results for the Tabula Muris FACS data set and generates the trait-cell type association plots for _seismic_ (Figure 3A), scDRS (Supplementary figure 6), FUMA (Supplementary figure 7), and S-MAGMA (Supplementary figure 8). 