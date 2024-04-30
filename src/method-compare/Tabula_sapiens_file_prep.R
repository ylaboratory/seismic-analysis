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
assay(ts_obj, "cpm") <- scuttle::calculateCPM(ts_obj, assay.type = "counts", size.factors = colSums(assay(ts_obj, "counts")))

#calculate mean expression and map to human gene id
ts_mean <- calc_ct_mean(ts_obj, assay_name = "cpm", ct_label_col = "cluster_name")

ts_mean_hsa <- seismicGWAS::translate_gene_ids(t(ts_mean), from = "hsa_ensembl")

#print magma files
print_magma_fuma_tbl(t(ts_mean_hsa), "MAGMA", main_table_path = here("data","expr","Tabula_sapiens","new_ts.top10_magma.txt"),
                     aux_table_path = here("data","expr","Tabula_sapiens","new_ts.magma.aux.txt"))


###### fuma data preparation ##### 
#logcpm
ts_obj <- scuttle::logNormCounts(ts_obj, assay.type = "cpm",size_factors=rep(1, ncol(ts_obj))) #no size factors

#calculate mean expression and map to human gene id (change assay)
ts_mean_log <- calc_ct_mean(ts_obj, assay_name = "logcounts", ct_label_col = "cluster_name")

ts_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(ts_mean_log), from = "hsa_ensembl")

#print fuma files
print_magma_fuma_tbl(t(ts_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Tabula_sapiens","new_ts.fuma.txt"),
                     aux_table_path = here("data","expr","Tabula_sapiens","new_ts.fuma.aux.txt"))


##scdrs 
#cells to keep (only cell types in the same analysis)
ts_obj_clean <- ts_obj[,ts_obj$cluster_name %in% colnames(ts_mean_log_hsa)]

facs_ad <- AnnData(X = assay(ts_obj_clean,"counts")%>% set_rownames(rownames(rowData(ts_obj_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(ts_obj_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(ts_obj_clean)) %>% set_rownames(.$symbol) )

facs_ad %>%  write_h5ad(here("data","expr","Tabula_sapiens","ts.clean.h5ad"))

