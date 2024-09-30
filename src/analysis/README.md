# seismic Analysis Scripts

The `src/analysis` directory contains scripts for performing seismic analysis on various datasets. These scripts identify trait-associated cell type enrichments and conduct influential gene analysis for selected trait-cell type pairs.

## Prerequisites

- Expression datasets should be prepared in advance (after processing with `src/data-cleaning`).
- GWAS summary statistics should be processed and transformed to gene-level MAGMA z-scores. (See details in the data preprocessing section of the main README)

## Scripts and Their Functions

1. `Saunders_analysis.R`
   - Analyzes brain single-cell RNA-seq data from Saunders et al. (2018)
   - Identifies cell types associated with neurodegenerative disorders and related traits
   - Performs analysis at multiple granularities

2. `Saunders_inf_analysis.R`
   - Conducts influential gene analysis on results from `Saunders_analysis.R`
   - Performs GO enrichment analysis on influential genes for selected trait-cell type pairs

3. `Tabula_muris_analysis.R`
   - Analyzes mouse single-cell RNA-seq data from Tabula Muris Consortium (2018)
   - Covers 20 organs to identify cell types associated with various traits
   - Separately analyzes FACS and droplet datasets

4. `Tabula_muris_ws.R`
   - Analyzes Tabula Muris FACS datasets using the same trait GWAS summary statistics
   - Employs varying window sizes as MAGMA gene analysis parameters

5. `Tabula_sapiens_analysis.R`
   - Analyzes human single-cell RNA-seq data from Tabula Sapiens Consortium (2019)