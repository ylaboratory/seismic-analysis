Several files are included in the `tools` directory to facilitate the analysis. These scripts are used for data processing, annotation, and analysis. The following is a brief description of each script:

- `annotate_snp.sh`: This script is used for annotating single nucleotide polymorphisms (SNPs) in genomic data. It takes SNP data as input and adds relevant annotations such as gene names, functional impact, and allele frequencies. 

- `fuma_magma_batch.sh`: This script is used for running batch jobs of the FUMA MAGMA tool. FUMA MAGMA is a gene-based association analysis tool that identifies genes associated with a trait or disease using summary statistics from genome-wide association studies (GWAS).

- `magma_fuma_file_prep.R`: This script prepares the input files for S-MAGMA (MAGMA cell-type-specific gene set enrichment analysis) and FUMA, which are further taken as the input for the MAGMA software. 

- `magma_gene_zscore_analysis.sh`: This script transforms GWAS summary statistics files to gene-level z-scores for gene-based association analysis using MAGMA. To ensure the correct usage, the MAGMA software should be downloaded and installed on the system. 

- `read_dropviz_data.R`: This script reads and processes data from the DropViz database. The orignal code did not functiona as expected and was modified to work with the data.

- `sparse_mat_util.R`: This script contains utility functions for working with sparse matrices, which are commonly used in single-cell RNA-seq data analysis. 