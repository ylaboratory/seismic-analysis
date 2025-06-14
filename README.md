# Analysis scripts for _seismic_

_seismic_ is a computationally efficient, powerful, and interpretable approach for
identifying associations between complex traits and cell type-specific expression
using GWAS summary statistics and single-cell RNA-seq data.
This repository contains all code and scripts required for reproducing the results
and figures of the paper,
[Disentangling associations between complex traits and cell types with _seismic_](https://www.biorxiv.org/content/10.1101/2024.05.04.592534).

All data for reproducing the analyses can be found on [zenodo](https://zenodo.org/records/15582078).
We recommend downloading these files and unzipping into this base directory before
proceeding. All file paths will assume that downloaded `all_data` file from zenodo is unzipped
and placed under this repo's parent directory.

Scripts are located in the `src` folder and are divided into the following directories:

- `analysis`: scripts for performing seismic analysis on various datasets.
- `causal-sim`: scripts for simulating trait-specific gene expression in cell groups and evaluating the performance of different cell type-trait association detection methods
- `data-cleaning`: scripts for preparing the expression data sets from raw data downloaded directly from their respective sources
- `extra`: scripts for extra experiments that analyze the performance of seismic
- `method-compare`: scripts for running other methods for baseline comparison 
- `null-sim`: scripts for performing null simulations to assess the calibration of different cell type-trait association frameworks
- `plot-figure`: scripts for generating all the figures in the paper
- `runtime`: scripts for evaluating the runtime performance of seismic and other cell type-trait association frameworks
- `tools`: accessory scripts to facilitate data processing, annotation, and analysis

## Basic setup

The _seismic_ framework is available as an R package
[_seismicGWAS_](https://github.com/ylaboratory/seismic).
Please refer to the package and accompanying vignette for more about usage,
requirements, and installation. To quickly install the package in R
use `devtools`:

```R
devtools::install_github("ylaboratory/seismicGWAS")
```

Other requirements including other tool prerequisites are located in the
README files in respective subdirectories.

We refer users to the [_seismicGWAS_](https://github.com/ylaboratory/seismic) repo
for detailed instructions on how to preprocess your own data
and run _seismic_ on your own data.
