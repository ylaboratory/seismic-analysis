# Analysis Tools
The `tools` directory contains various scripts to facilitate data processing, annotation, and analysis. These scripts are essential components of the overall analysis pipeline. Below is a brief description of each script:

- `annotate_snp.sh`: Annotates single nucleotide polymorphisms (SNPs) in genomic data.
  - Input: SNP data
  - Output: Annotated SNPs with gene names, functional impact, and allele frequencies

- `fuma_magma_batch.sh`: Runs batch jobs for the FUMA MAGMA tool.
  - Purpose: Performs gene-based association analysis using GWAS summary statistics
  - Identifies genes associated with traits or diseases

- `magma_fuma_file_prep.R`: Prepares input files for S-MAGMA and FUMA analyses.
  - Output: Files ready for input into MAGMA software
  - Used for cell-type-specific gene set enrichment analysis

- `magma_gene_zscore_analysis.sh`: Transforms GWAS summary statistics to gene-level z-scores.
  - Purpose: Enables gene-based association analysis using MAGMA
  - Note: Requires MAGMA software to be installed on the system

- `read_dropviz_data.R`: Reads and processes data from the DropViz database.
  - Note: Original code was modified to function correctly with the data

- `sparse_mat_util.R`: Contains utility functions for working with sparse matrices.
  - Application: Commonly used in single-cell RNA-seq data analysis