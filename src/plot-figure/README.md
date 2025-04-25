# Generate all figures in the paper using the result files
The `plot-figure` directory contains the scripts for generating all the figures in the paper. The scripts are organized as follows:

`null_sim_plot.R`
   - Processes null simulation results for scDRS, FUMA, and S-MAGMA
   - Generates calibration plots (Figure 2A, Supplementary Figure 1)

`causal_sim_plot.R`
   - Processes causal simulation results for all four frameworks
   - Generates power plots (Figure 2B, Supplementary Figures 2-3)
   - Produces influential gene precision-recall plots (Figure 2C, Supplementary Figure 4)

`runtime_plot.R`
   - Processes runtime evaluation results
   - Generates runtime comparison plots (Figure 2D)

`analysis_TM_FACS_plot.R`
   - Processes Tabula Muris FACS dataset analysis results
   - Generates:
     - Trait-cell type association plots for seismic (Figure 3A), scDRS (Supp. Fig. 6), FUMA (Supp. Fig. 7), and S-MAGMA (Supp. Fig. 8)
     - Difference heatmaps between seismic and other methods (Supp. Figs. 9-11)
     - Venn diagrams of significant associations (Figure 3B)
     - Cross-method correlation plots (Figure 3D, Supp. Fig. 12)

`varied_ws_plot.R`
   - Analyzes MAGMA gene-level Z-scores with different window sizes
   - Generates correlation plots (Supplementary Figure 5)

`analysis_TM_droplet_plot.R`
   - Processes Tabula Muris droplet dataset results
   - Generates:
     - Trait-cell type association plots for seismic (Supp. Fig. 13)
     - Cross-method correlation plots (Supp. Fig. 16)

`analysis_TS_plot.R`
   - Processes Tabula Sapiens dataset results
   - Generates:
     - Trait-cell type association plots for seismic (Supp. Fig. 14)
     - Cross-method correlation plots (Supp. Fig. 17)

`cross_ds_plot.R`
   - Maps cell type labels across datasets
   - Generates cross-dataset correlation plots (Supp. Fig. 15)

`analysis_Saunders_plot.R`
   - Processes Saunders et al. mouse brain scRNA-seq dataset results
   - Generates:
     - PD-associated cell types across granularities (Figure 4)
     - AD-associated cell types for different GWAS endpoints (Figure 5AB)

1`inf_analysis_plot.R`
    - Processes influential gene analysis results
    - Generates:
      - Influential gene plots for specific trait-cell type associations (Figure 5CD, left)
      - Enriched GO pathway plots (Figure 5CD, right)