# Generate all figures in the paper using the result files
The `plot-figure` directory contains the scripts for generating all the figures in the paper. The scripts are organized as follows:

- `null_sim_plot.R` processes the null simulation results for the scDRS, FUMA, and S-MAGMA frameworks and generates the calibration plots (Figure 2A, Supplementary Figure 1) for the three frameworks.

- `causal_sim_plot.R` processes the causal simulation results for the four frameworks in various conditions and generates:
    - The power plots (Figure 2B, Supplementary Figure 2-3). 
    - The influential gene precision recall plots (Figure 2C, Supplementary figure 4).

- `runtime_plot.R` processes the runtime evaluation results for the four frameworks and generates the runtime plots (Figure 2D) for the four frameworks. 

- `analysis_TM_FACS_plot.R` processes the analysis results for the Tabula Muris FACS dataset and generates the fllowing plots: 
    - Trait-cell type association plots for _seismic_ (Figure 3A), scDRS (Supplementary figure 6), FUMA (Supplementary figure 7), and S-MAGMA (Supplementary figure 8).
    - Difference heatmaps for all trait-cell type pairs between _seismic_ and any alternative method (Supplementary figure 9-11).
    - Venn diagrams for the number of significant association pairs for _seismic_, scDRS, FUMA, and S-MAGMA (Figure 3B).
    - Cross-method correlation plots for all traits (Figure 3D, Supplementary figure 12).

- `varied_ws_plot.R` collecting cell-type association using MAGMA gene-level Z-score with different window size arguments and generates correlation between the results (Supplementary figure 5).

- `analysis_TM_droplet_plot.R` processes the analysis results for the Tabula Muris droplet dataset and generates the fllowing plots: 
    - Trait-cell type association plots for _seismic_ (Supplementary figure 13).
    - Cross-method correlation plots for all traits (=Supplementary figure 16).

- `analysis_TS_plot.R` processes the analysis results for the Tabula sapiens dataset and generates the fllowing plots: 
    - Trait-cell type association plots for _seismic_ (Supplementary figure 14).
    - Cross-method correlation plots for all traits (Supplementary figure 17).

- `cross_ds_plot.R` maps the cell type labels between different datasets (while kepping cell types with only similar tissues and cell ontology) and generates the cross-dataset correlation plots (Supplementary figure 15).

- `analysis_Saunders_plot.R` processes the analysis results for the Saunders et al. mouse whole brain scRNA-seq dataset and generates the fllowing plots: 
    - Parkinson's disease(PD) associated cell types across multiple analysis granularities using _seismic_ (Figure 4).
    - Plots for top associated cell type for Alzheimer's disease(AD) GWAS with different endpoints (clinical AD diagnosis and cerebrospinal fluid tau level) using _seismic_ (Figure 5AB).

