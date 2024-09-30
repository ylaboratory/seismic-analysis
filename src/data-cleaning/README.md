# Data cleaning scripts

This folder contains scripts for preparing the expression data sets from raw data downloaded directly from the respective sources. The scripts are organized as follows:

- `Saunders_data_cleaning.R` processes the [raw data from Saunders et al. (2018)](ttp://dropviz.org/). It filters out cells labeled as singleton, doublet, or outliers. Subcluster annotations are manually curated to remove typos and inconsistencies, which are used for the fine_cluster granularity. Cluster labels are extracted from the supplementary file of the original paper. Neuron subclass labels are annotated based on the primary neurotransmitter type.

- `Tabula_muris_data_cleaning.R` cleans the [raw data from Tabula Muris Consortium (2018)](https://figshare.com/projects/Tabula_Muris_Transcriptomic_characterization_of_20_organs_and_tissues_from_Mus_musculus_at_single_cell_resolution/27733.). The script downloads data from the figshare repository and extracts expression data from the Seurat object.

- `Tabula_sapiens_data_cleaning.R` processes the [raw data from Tabula Sapiens Consortium (2022)]( https://figshare.com/projects/Tabula_Sapiens/100973). It extracts expression data from the Seurat object and maps cell ontology classes to cell ontology IDs using a separate Python script.

- `Tabula_sapiens_ontology_cleaning.py` maps cell ontology classes to cell ontology IDs for the Tabula Sapiens dataset. It uses [the whole cell ontology mapping file](https://obofoundry.org/ontology/cl.html). Any annotation that cannot be mapped is manually curated.

- `gwas_processing.R` processes GWAS summary statistics to MAGMA gene-level z-scores. It downloads GWAS summary statistics from sources listed in the supplementary file, calculates MAGMA gene-level z-scores, and saves .tsv files in the `data/gwas` directory for scDRS analysis.

- `sim_data_processing.R` generates random seeds for null and causal simulations. It creates shuffled expression data, gene sets, and random cell samples. It also samples 10 random traits to generate random MAGMA gene-level z-scores for null simulation. Output is saved in `data/expr` and `data/sim` directories.

- `runtime_data_processing.R` prepares data for runtime evaluation. It samples expression datasets of different sizes, schedules tasks using podman, measures execution time for each framework, and saves runtime data in the `data/runtime` directory.
