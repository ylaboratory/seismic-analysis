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
assay(facs_obj_sce, "cpm")  <- scuttle::calculateCPM(facs_obj_sce, assay.type = "counts", size.factors = colSums(assay(facs_obj_sce, "counts")))
assay(droplet_obj_sce, "cpm")  <- scuttle::calculateCPM(droplet_obj_sce, assay.type = "counts", size.factors = colSums(assay(droplet_obj_sce, "counts")))

#calculate mean expression and map to human gene id
facs_mean <- calc_ct_mean(facs_obj_sce, assay_name = "cpm", ct_label_col = "cluster_name")
droplet_mean <- calc_ct_mean(droplet_obj_sce, assay_name = "cpm", ct_label_col = "cluster_name")

facs_mean_hsa <- seismicGWAS::translate_gene_ids(t(facs_mean), from = "mmu_symbol")
droplet_mean_hsa <- seismicGWAS::translate_gene_ids(t(droplet_mean), from = "mmu_symbol")

#print magma files
print_magma_fuma_tbl(t(facs_mean_hsa), "MAGMA", main_table_path = here("data","expr","Tabula_muris","new_tm_facs.top10_magma.txt"),
                     aux_table_path = here("data","expr","Tabula_muris","new_tm_facs.magma.aux.txt"))
print_magma_fuma_tbl(t(droplet_mean_hsa), "MAGMA", main_table_path = here("data","expr","Tabula_muris","new_tm_droplet.top10_magma.txt"),
                     aux_table_path = here("data","expr","Tabula_muris","new_tm_droplet.magma.aux.txt"))


###### fuma data preparation ##### 
#logcpm
facs_obj_sce <- scuttle::logNormCounts(facs_obj_sce, assay.type = "cpm", size_factors=rep(1, ncol(facs_obj_sce))) #no size factors
droplet_obj_sce <- scuttle::logNormCounts(droplet_obj_sce, assay.type = "cpm",size_factors=rep(1, ncol(droplet_obj_sce))) #no size factors

#calculate mean expression and map to human gene id (change assay)
facs_mean_log <- calc_ct_mean(facs_obj_sce, assay_name = "logcounts", ct_label_col = "cluster_name")
droplet_mean_log <- calc_ct_mean(droplet_obj_sce, assay_name = "logcounts", ct_label_col = "cluster_name")

facs_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(facs_mean_log), from = "mmu_symbol")
droplet_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(droplet_mean_log), from = "mmu_symbol")

#print fuma files
print_magma_fuma_tbl(t(facs_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Tabula_muris","new_tm_facs.fuma.txt"),
                     aux_table_path = here("data","expr","Tabula_muris","new_tm_facs.fuma.aux.txt"))
print_magma_fuma_tbl(t(droplet_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Tabula_muris","new_tm_droplet.fuma.txt"),
                     aux_table_path = here("data","expr","Tabula_muris","new_tm_droplet.fuma.aux.txt"))

##scdrs 
#cells to keep (only cell types in the same analysis)
facs_obj_sce_clean <- facs_obj_sce[,facs_obj_sce$cluster_name %in% colnames(facs_mean_hsa)]
droplet_obj_sce_clean <- droplet_obj_sce[,droplet_obj_sce$cluster_name %in% colnames(droplet_mean_hsa)]

facs_ad = AnnData(X = assay(facs_obj_sce_clean,"counts")%>% set_rownames(rownames(rowData(facs_obj_sce_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(facs_obj_sce_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(facs_obj_sce_clean)) %>% set_rownames(.$symbol) )
droplet_ad = AnnData(X = assay(droplet_obj_sce_clean,"counts")%>% set_rownames(rownames(rowData(droplet_obj_sce_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(droplet_obj_sce_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(droplet_obj_sce_clean)) %>% set_rownames(.$symbol) )

facs_ad %>%  write_h5ad(here("data","expr","Tabula_muris","facs.clean.h5ad"))
droplet_ad %>%  write_h5ad(here("data","expr","Tabula_muris","droplet.clean.h5ad"))



