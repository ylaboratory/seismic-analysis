# _seismic_ Analysis Scripts

This directory contains scripts for performing _seismic_ analyses
on various datasets. These scripts identify trait-associated cell
type enrichments and conduct influential gene analysis for selected
trait-cell type pairs.

## Optional Prerequisites

We have currently written most scripts to use data in the
[seismic data repo](https://zenodo.org/records/15582078).

If building from the raw data sources you must make sure that:

- Expression datasets are prepared in advance (after processing with `src/data-cleaning`)
- GWAS summary statistics are processed and transformed to gene-level MAGMA z-scores

## Script Index

1. `Saunders_analysis.R`
   - Analyzes brain single-cell RNA-seq data from Saunders et al. (2018)
   - Identifies cell types associated with neurodegenerative disorders and related traits
   - Performs analysis at multiple granularities
   - Conducts influential gene analysis
   - Performs GO enrichment analysis on influential genes for selected trait-cell type pairs

2. `Tabula_muris_analysis.R`
   - Analyzes mouse single-cell RNA-seq data from Tabula Muris Consortium (2018)
   - Covers 20 organs to identify cell types associated with various traits
   - Separately analyzes FACS and droplet datasets

3. `Tabula_muris_ws.R`
   - Analyzes Tabula Muris FACS datasets using the same trait GWAS summary statistics
   - Employs varying window sizes as MAGMA gene analysis parameters
   - _Note: MAGMA files of varying window sizes need to be first generated with `tools/magma_gene_zscore_analysis.sh`_

4. `Tabula_sapiens_analysis.R`
   - Analyzes human single-cell RNA-seq data from Tabula Sapiens Consortium (2019)