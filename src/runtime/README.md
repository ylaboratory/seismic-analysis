# Runtime Evaluation
The `runtime` subdirectory contains the code for evaluating the runtime performance of _seismic_ analysis procedures. It includes scripts for measuring the execution time for the four frameworks on the  scRNA-seq datasets with differnet sizes and the script to schedule all the tasks using podman.

- `podman_run_per_ds.sh` schedules the runtime evaluation tasks for the four frameworks sequentially on the same scRNA-seq dataset processed with `src/data-cleaning/sim_data_processing.R` and save the results in the `results/runtime` directory. For each of the four frameworks, the script will call the corresponding script or command to run the analysis and measure the execution time. The procedures are:

    - `scDRS`: By calling the `scdrs_runtime.py` script, which takes the expression dataset and several gene sets as input, the script will record the runtime for preprocessing the dataset, calculating the scDRS scores and evaluating cell-type level trait associations. The output from different files will be merged with the `munge_runtime_scdrs.py` script.
    - `FUMA`: `fuma_data_gene.R` generates the input files of the current expression dataset for the FUMA model and further calls the MAGMA software to calculate the cell-type level trait associations. The runtime is recorded by the script.
    - `S-MAGMA`: Similar to the FUMA procedure, the script will call the `magma_data_gene.R` script to generate the input files for the S-MAGMA model and record the runtime for the analysis.
    - `seismic`: The `seismic_runtime.py` script will take the expression dataset and the gene sets as input and record the runtime for calculating the specificity score and evaluating cell-type level trait associations.
