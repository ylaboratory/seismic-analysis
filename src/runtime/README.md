# Runtime Evaluation

This directory contains scripts for evaluating the runtime performance
of _seismic_ and other cell type-trait association frameworks. The analysis
is performed on scRNA-seq datasets of varying sizes, processed
using `src/data-cleaning/sim_data_processing.R`.

The main script and its components are as follows:

- `podman_run_per_ds.sh`: A shell script that schedules runtime evaluation tasks for all frameworks sequentially on each scRNA-seq dataset. It calls framework-specific scripts and commands, measures execution time, and saves results in the `results/runtime` directory. The evaluation process for each framework is as follows:
- `scdrs_runtime.py`:  measures preprocessing time, scDRS score calculation time, cell-type level trait association evaluation time
- `munge_runtime_scdrs.py`: merges the scDRS results
- `fuma_data_gene.R`: generates input files for the FUMA model, calls MAGMA software for cell-type level trait associations
- `magma_data_gene.R`: similar to FUMA, generates input files for the S-MAGMA model
- `seismic_runtime.py`: measures time for specificity score calculation and cell-type level trait association evaluation

All runtime results are saved in the `results/runtime` directory for further analysis and comparison.