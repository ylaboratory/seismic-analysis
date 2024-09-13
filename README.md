# Analysis scripts for *seismic*
*seismic* is is a computationally efficient, powerful, and interpretable approach for identifying associations between complex traits and cell type-specific expression using GWAS summary statistics and single-cell RNA-seq data. This repository contains all code and scripts required for reproducing the results and figures of the paper, [Disentangling associations between complex traits and cell types with *seismic*](https://www.biorxiv.org/content/10.1101/2024.05.04.592534).  The repository is organized as follows:

- bin: This directory stores applications, such as MAGMA and scDRS.

- data: All processed expression and GWAS data will be put here, as well as the intermediate files.

- results: The results of all the analyses, including all real dataset analyses, simulation, runtime analysis, and benchmarking.

- raw: This directory contains the raw and unprocessed data.

- ref: Reference files and datasets that are used as inputs for the analysis are stored in this directory. Such as MAGMA auxillary files. 

- src: Source code files for the analysis, including all scripts for data generation, data analysis, tool scripts and figure generation.

## Environment set up
The *seismic* framework is packed up as an R package named [*seismicGWAS*](https://github.com/ylaboratory/seismic) that is avialble for installation. Please refer to the package link and [vignette](https://ylaboratory.github.io/seismic/seismicGWAS.html) to know more about the usage and the download of the package. To download and install the package:

```{r}
devtools::install_github("ylaboratory/seismicGWAS")
```
This will be sufficient for the *seismic* trait-associated cell type analysis and influential gene analysis. To run the benchmarking analysis, please refer to the script in the `src/method-compare` subdirectory.
