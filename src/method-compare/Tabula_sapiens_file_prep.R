#prepare files for scDRS/FUMA/MAGMA 

##### 1. load packages and data#######
###load packages
if (!require("here")){
  install.packages("here")
  library("here")
}

if (!require("magrittr")){
  install.packages("magrittr")
  library("magrittr")
}
if (!require("tidyverse")){
  install.packages("tidyverse")
  library("tidyverse")
}

if (!require("SingleCellExperiment")){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("SingleCellExperiment")
  library("SingleCellExperiment")
} 

if (!require("seismicGWAS")){
  if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

#create anndata object and print out
if (!require("anndata")){
  install.packages("anndata")
  library("anndata")
}

#load function
source(here("src","tools","magma_fuma_file_prep.R"))
source(here("src","tools","sparse_mat_util.R"))

#load the object that's filtered already in the Tabula_sapeins_analysis scripts
load(here("data","expr","Tabula_sapiens","TS_processed.rda"))

###### magma data preparation #####
#cpm
assay(ts_obj, "cpm")  = scuttle::calculateCPM(ts_obj, assay.type = "counts")

#calculate mean expression and map to human gene id
data("mmu_hsa_mapping")
ts_obj = cal_stat(ts_obj, meta_data = as.data.frame(colData(ts_obj)), group = "cluster_name",assay_name ="cpm" ,mean_only=T)
ts_obj= trans_mmu_to_hsa_stat(ts_obj, gene_mapping_table = mmu_hsa_mapping, from = "hsa_ensembl", to = "hsa_entrez")

#print magma files
print_magma_fuma_tbl(ts_obj, "MAGMA", main_table_path = here("data","expr","Tabula_sapiens","ts.top10_magma.txt"),
                     aux_table_path = here("data","expr","Tabula_sapiens","ts.magma.aux.txt"))

###### fuma data preparation ##### 
#logcpm
ts_obj = scuttle::logNormCounts(ts_obj, assay.type = "cpm",size_factors=rep(1, ncol(ts_obj))) #no size factors

#calculate mean expression and map to human gene id (change assay)
ts_obj = cal_stat(ts_obj, meta_data = as.data.frame(colData(ts_obj)), group = "cluster_name",assay_name ="logcounts", mean_only=T)
ts_obj= trans_mmu_to_hsa_stat(ts_obj, gene_mapping_table = mmu_hsa_mapping, from = "hsa_ensembl", to = "hsa_entrez")

#print fuma files
print_magma_fuma_tbl(ts_obj, "FUMA", main_table_path = here("data","expr","Tabula_sapiens","ts.fuma.txt"),
                     aux_table_path = here("data","expr","Tabula_sapiens","ts.fuma.aux.txt"))

##scdrs 
#cells to keep (only cell types in the same analysis)
ts_obj_clean = ts_obj[,ts_obj$cluster_name %in% names(get_seismic_ct_info(ts_obj, "cell_num"))]

facs_ad = AnnData(X = assay(ts_obj_clean,"counts")%>% set_rownames(rownames(rowData(ts_obj_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(ts_obj_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(ts_obj_clean)) %>% set_rownames(.$symbol) )

facs_ad %>%  write_h5ad(here("data","expr","Tabula_sapiens","ts.clean.h5ad"))

