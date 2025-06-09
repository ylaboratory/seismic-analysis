# Runtime Evaluation

This directory contains scripts for evaluating the runtime performance of seismic and other cell type-trait association frameworks. The analysis is performed on scRNA-seq datasets of varying sizes, processed using `src/data-cleaning/sim_data_processing.R`.

The main script and its components are as follows:

- `podman_run_per_ds.sh`: A shell script that schedules runtime evaluation tasks for all frameworks sequentially on each scRNA-seq dataset. It calls framework-specific scripts and commands, measures execution time, and saves results in the `results/runtime` directory. The evaluation process for each framework is as follows:

- scDRS:
  - Script: `scdrs_runtime.py`
  - Input: Expression dataset and gene sets
  - Measures: Preprocessing time, scDRS score calculation time, cell-type level trait association evaluation time
  - Post-processing: Results merged using `munge_runtime_scdrs.py`

- FUMA:
  - Script: `fuma_data_gene.R`
  - Process: Generates input files for the FUMA model, calls MAGMA software for cell-type level trait associations
  - Measures: Total runtime for the entire process

- S-MAGMA:
  - Script: `magma_data_gene.R`
  - Process: Similar to FUMA, generates input files for the S-MAGMA model
  - Measures: Total runtime for the analysis

- _seismic_:
  - Script: `seismic_runtime.py`
  - Input: Expression dataset and gene sets
  - Measures: Time for specificity score calculation and cell-type level trait association evaluation

All runtime results are saved in the `results/runtime` directory for further analysis and comparison.