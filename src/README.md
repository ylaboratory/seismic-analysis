# Source code

This directory includes most of the source code for the projects (e.g., data analysis, data cleaning, benchmarking, null simulation, runtime evaluation, and extra tool functions). Here is the organization of the code:

### Data Analysis
The `analysis` subdirectory contains the main code for analyzing the scRNA-seq data with *seismic*. The data are previously cleaned with scripts in the `data-cleaning` subdirectory. The analysis code includes scripts for normalization, quality control, *seismic* cell-type level trait association analysis. The code are organized based on the scRNA-seq data sets. 

### Data Cleaning
The `data-cleaning` subdirectory contains the code for cleaning and preprocessing the scRNA-seq data from the `raw` directory. All labels are refined and the cells with confusing labels are removed (such as hepatocytes in heart for the TS data set). The data are then saved in the `data/expr` directory for further analysis. Scripts for generating data sets (scRNA-seq, GWAS, and random seeds) for the null simulation and runtime analysis are also included in this subdirectory.

### Method Comparison
The `method-comparison` subdirectory contains the code for generating benchmarking results for different frameworks (scDRS, FUMA, S-MAGMA). It includes data preparation steps whcih output the required file format for these frameworks and scripts for running these framework pipelines. The results are saved in the `results` directory for comparison. 

### Null Simulation
The `null-simulation` subdirectory contains the code for assessing each framework's calibration for the null hypothesis. It includes scripts for running the analysis pipelines and reporting the significance (P-values) across 10,000 samples of random assembled target cell types. The shuffled expression and gene sets are previously prepared in the `data-cleaning` subdirectory.

### Causal Simulation
The `causal-simulation` subdirectory contains the code for assessing the causal inference performance of the *seismic* framework. It includes scripts for running the analysis pipelines and reporting the statistical significance (or plus influential gene analysis for *seismic*). The synthetic data are generated with the `data-cleaning/causal_sim_data_scDesign3.R` script.

### Runtime Evaluation
The `runtime-evaluation` subdirectory contains the code for evaluating the runtime performance of the seismic analysis algorithms. It includes scripts for measuring the execution time for the four frameworks on the same scRNA-seq data set and the script to schedule all the tasks using podman. 

### Extra Experiment and Plot
The `extra` subdirectory contains the code for generating extra figures and results that are not included in the main text. It includes scripts for analyzing the false discovery rate of the *seismic* framework, label perturbation analysis, and other supplementary figures. The results are saved in the `results` directory.

### Extra Tool Functions
The `tools` subdirectory contains additional utility functions and tools that are used across different parts of the *seismic* projects. 