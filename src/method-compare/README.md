# Analyzing the data with the other frameworks

- `run_all_frameworks.sh` is the shell scripts containing the commands to run the scDRS, FUMA, and S-MAGMA frameworks on the same scRNA-seq data set.
The script is used to generate the benchmarking results for the method comparison.
The results are saved in the `results` directory for comparison. Before running these commands, the data should be processed into the required file formats, as described in the other scripts in the `method-comparison` subdirectory.

- `Tabula_muris_file_prep.R` is the R script for processing the expression data for the Tabula Muris dataset (including both the FACS and droplet datasets).
The script reads the expression data saved after analysis as in `src/analysis/Tabula_muris_analysis.R` and processes it into the required file format for the scDRS, FUMA, and S-MAGMA frameworks.
The processed data is saved in the `data/expr/Tabula_muris` directory for the method comparison.

- `Tabula_sapiens_file_prep.R` is similar to `Tabula_muris_file_prep.R` but for the Tabula Sapiens dataset. T The processed data is saved in the `data/expr/Tabula_sapiens`.

- `Saunders_file_prep.R` is similar to `Tabula_muris_file_prep.R` but for the Saunders dataset. The data is processed using 5 different analysis granularity and is saved in the `data/expr/Saunders` directory.

