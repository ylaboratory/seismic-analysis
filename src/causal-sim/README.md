# Causal Simulation for Cell Type-Trait Association Methods

This directory contains scripts for simulating trait-specific gene
expression in cell groups and evaluating the performance of different
cell type-trait association detection methods.
We use the [scDesign3](https://www.nature.com/articles/s41587-023-01772-1)
approach to generate synthetic data.

## Workflow

1. `sample_cells.R`: Sample cells from datasets to create cell groups for perturbation
2. `fit_scDesign3.R`: Fit scDesign3 models for datasets with new cell type information
3. `merge_scDesign3_model_obj.R`: merge individual scDesign3 models into a single object
4. `causal_sim_data_gene.R`: simulate expression profiles
5. `causal_sim.R`: run _seismic_, S-MAGMA, and FUMA and assess performance
6. `scdsrs_causal_sim.py`: run scDRS and assess performance

## Full Pipeline Execution

The `causal_sim_all.sh` shell script contains all commands to:
- Fit scDesign3 models
- Run `causal_sim.R` and `scdrs_causal_sim.py`
- Assess performance

_Note: `causal_sim.R` and `scdsrs_causal_sim.py` are designed to run in parallel for faster processing._
