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

#load the object that's filtered already in the Tabula_muris_analysis scripts
load(here("data","expr","Tabula_muris","TM_processed.rda"))

###### magma data preparation #####
#cpm
assay(facs_obj_sce, "cpm")  = scuttle::calculateCPM(facs_obj_sce, assay.type = "counts")
assay(droplet_obj_sce, "cpm")  = scuttle::calculateCPM(droplet_obj_sce, assay.type = "counts")

#calculate mean expression and map to human gene id
facs_obj_sce = cal_stat(facs_obj_sce, meta_data = as.data.frame(colData(facs_obj_sce)), group = "cluster_name",assay_name ="cpm" ,mean_only=T)
droplet_obj_sce = cal_stat(droplet_obj_sce, meta_data = as.data.frame(colData(droplet_obj_sce)), group = "cluster_name",assay_name ="cpm" ,mean_only=T)

facs_obj_sce= trans_mmu_to_hsa_stat(facs_obj_sce, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
droplet_obj_sce= trans_mmu_to_hsa_stat(droplet_obj_sce, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")

#print magma files
print_magma_fuma_tbl(facs_obj_sce, "MAGMA", main_table_path = here("data","expr","Tabula_muris","tm_facs.top10_magma.txt"),aux_table_path = here("data","expr","Tabula_muris","tm_facs.magma.aux.txt"))
print_magma_fuma_tbl(droplet_obj_sce, "MAGMA", main_table_path = here("data","expr","Tabula_muris","tm_droplet.top10_magma.txt"),aux_table_path = here("data","expr","Tabula_muris","tm_droplet.magma.aux.txt"))


###### fuma data preparation ##### 
#logcpm
facs_obj_sce = scuttle::logNormCounts(facs_obj_sce, assay.type = "cpm",size_factors=rep(1, ncol(facs_obj_sce))) #no size factors
droplet_obj_sce = scuttle::logNormCounts(droplet_obj_sce, assay.type = "cpm",size_factors=rep(1, ncol(droplet_obj_sce))) #no size factors

#calculate mean expression and map to human gene id (change assay)
facs_obj_sce = cal_stat(facs_obj_sce, meta_data = as.data.frame(colData(facs_obj_sce)), group = "cluster_name",assay_name ="logcounts", mean_only=T)
droplet_obj_sce = cal_stat(droplet_obj_sce, meta_data = as.data.frame(colData(droplet_obj_sce)), group = "cluster_name",assay_name ="logcounts", mean_only=T)

facs_obj_sce= trans_mmu_to_hsa_stat(facs_obj_sce, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
droplet_obj_sce= trans_mmu_to_hsa_stat(droplet_obj_sce, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")

#print fuma files
print_magma_fuma_tbl(facs_obj_sce, "FUMA", main_table_path = here("data","expr","Tabula_muris","tm_facs.fuma.txt"),aux_table_path = here("data","expr","Tabula_muris","tm_facs.fuma.aux.txt"))
print_magma_fuma_tbl(droplet_obj_sce, "FUMA", main_table_path = here("data","expr","Tabula_muris","tm_droplet.fuma.txt"),aux_table_path = here("data","expr","Tabula_muris","tm_droplet.fuma.aux.txt"))

##scdrs 
#cells to keep (only cell types in the same analysis)
facs_obj_sce_clean = facs_obj_sce[,facs_obj_sce$cluster_name %in% names(get_seismic_ct_info(facs_obj_sce, "cell_num"))]
droplet_obj_sce_clean = droplet_obj_sce[,droplet_obj_sce$cluster_name %in% names(get_seismic_ct_info(droplet_obj_sce, "cell_num"))]

facs_ad = AnnData(X = assay(facs_obj_sce_clean,"counts")%>% set_rownames(rownames(rowData(facs_obj_sce_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(facs_obj_sce_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(facs_obj_sce_clean)) %>% set_rownames(.$symbol) )
droplet_ad = AnnData(X = assay(droplet_obj_sce_clean,"counts")%>% set_rownames(rownames(rowData(droplet_obj_sce_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(droplet_obj_sce_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(droplet_obj_sce_clean)) %>% set_rownames(.$symbol) )

facs_ad %>%  write_h5ad(here("data","expr","Tabula_muris","facs.clean.h5ad"))
droplet_ad %>%  write_h5ad(here("data","expr","Tabula_muris","droplet.clean.h5ad"))



