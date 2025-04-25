# Causal Simulation for Cell Type-Trait Association Methods

This directory (`src/causal-sim`) contains scripts for simulating trait-specific gene expression in cell groups and evaluating the performance of different cell type-trait association detection methods. We use the [scDesign3](https://www.nature.com/articles/s41587-023-01772-1) approach to generate synthetic data.

## Workflow

- **Sample Cells**
   - Script: `sample_cells.R`
   - Purpose: Sample cells from datasets to create cell groups for perturbation
   - Output: New cell type information (printed and saved)

- **Fit scDesign3 Models**
   - Script: `fit_scDesign3.R`
   - Usage: `Rscript fit_scDesign3.R <dataset_file_path> <cell_type_info_file_path> <output_directory>`
   - Purpose: Fit scDesign3 models for datasets with new cell type information
   - Output: 10 individual scDesign3 models

- **Merge scDesign3 Models**
   - Script: `merge_scDesign3_model_obj.R`
   - Purpose: Merge individual scDesign3 models into a single object

- **Simulate Expression Profiles**
   - Script: `causal_sim_data_gene.R`
   - Purpose: 
     - Sample target gene sets for perturbation
     - Reload saved scDesign3 models
     - Simulate new expression profiles
   - Output: 
     - Synthetic data matrices
     - Target gene set information

- **Run Methods and Assess Performance**
   - Scripts: 
     - `causal_sim.R` (for _seismic_, S-MAGMA, and FUMA)
     - `scdsrs_causal_sim.py` (for scDSRS)
   - Purpose: Run methods on simulated data and evaluate performance
   - Output: 
     - Results for all cell types
     - _seismic_'s influential gene analysis for the target cell type

## Full Pipeline Execution

The `causal_sim.sh` shell script contains all commands to:
- Fit scDesign3 models
- Run `causal_sim.R` and `scdrs_causal_sim.py`
- Assess performance

Note: `causal_sim.R` and `scdsrs_causal_sim.py` are designed to run in parallel for faster processing.
