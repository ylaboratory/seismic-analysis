# Source code

This directory includes most of the source code for the projects (e.g., data analysis, data cleaning, benchmarking, null simulation, runtime evaluation, and extra tool functions). Here is the organization of the code:

### Data Analysis
The `data-analysis` subdirectory contains the main code for analyzing the scRNA-seq data with *seismic*. The data are previously cleaned with scripts in the `data-cleaning` subdirectory. The analysis code includes scripts for normalization, quality control, *seismic* cell-type level trait association analysis. The code are organized based on the scRNA-seq data sets. 

### Data Cleaning
The `data-cleaning` subdirectory contains the code for cleaning and preprocessing the scRNA-seq data from the `raw` directory. All labels are refined and the cells with confusing labels are removed (such as hepatocytes in heart for the TS data set). The data are then saved in the `data/expr` directory for further analysis. Scripts for generating data sets (scRNA-seq, GWAS, and random seeds) for the null simulation and runtime analysis are also included in this subdirectory.

### Method Comparison
The `method-comparison` subdirectory contains the code for generating benchmarking results for different frameworks (scDRS, FUMA, S-MAGMA). It includes data preparation steps whcih output the required file format for these frameworks and scripts for running these framework pipelines. The results are saved in the `results` directory for comparison. 

### Null Simulation
The `null-simulation` subdirectory contains the code for assessing each framework's calibration for the null hypothesis. It includes scripts for running the analysis pipelines and reporting the significance (P-values) across 10,000 samples of random assembled target cell types. The shuffled expression and gene sets are previously prepared in the `data-cleaning` subdirectory.

### Runtime Evaluation
The `runtime-evaluation` subdirectory contains the code for evaluating the runtime performance of the seismic analysis algorithms. It includes scripts for measuring the execution time for the four frameworks on the same scRNA-seq data set and the script to schedule all the tasks using podman. 

### Extra Tool Functions
The `tools` subdirectory contains additional utility functions and tools that are used across different parts of the *seismic* projects. 

- `annotate_snp.sh`: This script is used for annotating single nucleotide polymorphisms (SNPs) in genomic data. It takes SNP data as input and adds relevant annotations such as gene names, functional impact, and allele frequencies.

- `fuma_magma_batch.sh`: This script is used for running batch jobs of the FUMA MAGMA tool. FUMA MAGMA is a gene-based association analysis tool that identifies genes associated with a trait or disease using summary statistics from genome-wide association studies (GWAS).

- `magma_fuma_file_prep.R`: This script prepares the input files for S-MAGMA (MAGMA cell-type-specific gene set enrichment analysis) and FUMA, which are further taken as the input for the MAGMA software. 

- `magma_gene_zscore_analysis.sh`: This script transforms GWAS summary statistics files to gene-level z-scores for gene-based association analysis using MAGMA. To ensure the correct usage, the MAGMA software should be downloaded and installed on the system. 

- `read_dropviz_data.R`: This script reads and processes data from the DropViz database. The orignal code did not functiona as expected and was modified to work with the data.

- `sparse_mat_util.R`: This script contains utility functions for working with sparse matrices, which are commonly used in single-cell RNA-seq data analysis. 