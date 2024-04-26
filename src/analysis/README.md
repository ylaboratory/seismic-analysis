The directory `src/analysis` contains the code for the *seismic* analysis of the data, including trait-associated cell type enrichment analysis and influential gene analysis for interested trait-cell type pairs. The expression datasets has been previously prepared and the GWAS summary statistics has been processed and transformed to gene-level MAGMA z-score (to see the detail). The analysis code is organized as follows:

- `Saunders_analysis.R` Analysis of the brain single-cell RNA-seq data from Saunders et al. (2018) to identify cell types associated with several neurodegerative disorders (and related traits). Multiple analysis granularities are used.

- `Saunders_inf_analysis.R ` Influential gene analysis and influential gene GO enrichment for several intersted trait-cell type pairs resulted from the `Saunders_analysis.R` script. 

-  `Tabula_muris_analysis.R` Analysis of the mouse single-cell RNA-seq data from Tabula Muris Consortium (2018), which are comprised of 20 organs, to identify cell types associated with multiple various traits. Analysis of the FACS and droplet datasets are performed separately in the script. 

- `Tabula_muris_ws.R` Analysis of the TM FACS data sets using the same trait GWAS summary statistics but with a set of varying window sizes as MAGMA gene analysis parameters.

- `Tabula_sapiens_analysis.R` Analysis of the human single-cell RNA-seq data from Tabula Sapiens Consortium (2019). 

