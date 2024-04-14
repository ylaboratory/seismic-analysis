# Analysis scripts for *seismic*

This repository contains all code and scripts required for the analysis for the seismic paper. Here is the organization of the repository:

- bin: This directory stores applications, such as MAGMA and scDRS.

- data: All processed expression and GWAS data will be put here, as well as the intermediate files.

- results: The results of all the analyses, including all real dataset analyses, simulation, runtime analysis, and benchmarking.

- raw: This directory contains the raw and unprocessed data.

- ref: Reference files and datasets that are used as inputs for the analysis are stored in this directory. Such as MAGMA auxillary files. 

- src: Source code files for the analysis, including all scripts for data generation, data analysis, tool scripts and figure generation.

## Environment set up
For the basic *seismic* analysis, only the package should be download and installed. Please refer to the package link and vignette to know more about the usage
and the download of the package. To download and install the *seismic* package:

```{r}
devtools::install_github("ylaboratory/seismicGWAS")
```
