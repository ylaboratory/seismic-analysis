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

#load the object that's filtered already in the Tabula_muris_analysis scripts
load(here("data","expr","Saunders","Saunders_processed_multi.rda"))

###### magma data preparation #####
#cpm
assay(brain_sce, "cpm")  = scuttle::calculateCPM(brain_sce, assay.type = "counts")
assay(brain_sce_subclass, "cpm") = scuttle::calculateCPM(brain_sce_subclass, assay.type = "counts")
assay(brain_sce_region_class, "cpm") = scuttle::calculateCPM(brain_sce_region_class, assay.type = "counts")
assay(brain_sce_region_subclass, "cpm") = scuttle::calculateCPM(brain_sce_region_subclass, assay.type = "counts")
assay(brain_sce_region_cluster, "cpm") = scuttle::calculateCPM(brain_sce_region_cluster, assay.type = "counts")

#calculate mean expression and map to human gene id
brain_sce = cal_stat(brain_sce, meta_data = as.data.frame(colData(brain_sce)), group = "fine_cluster",assay_name ="cpm" ,mean_only=T)
brain_sce_subclass = cal_stat(brain_sce_subclass, meta_data = as.data.frame(colData(brain_sce_subclass)), group = "subclass",assay_name ="cpm" ,mean_only=T)
brain_sce_region_class = cal_stat(brain_sce_region_class, meta_data = as.data.frame(colData(brain_sce_region_class)), group = "region_class",assay_name ="cpm" ,mean_only=T)
brain_sce_region_subclass = cal_stat(brain_sce_region_subclass, meta_data = as.data.frame(colData(brain_sce_region_subclass)), group = "region_subclass",assay_name ="cpm" ,mean_only=T)
brain_sce_region_cluster = cal_stat(brain_sce_region_cluster, meta_data = as.data.frame(colData(brain_sce_region_cluster)), group = "region_cluster",assay_name ="cpm" ,mean_only=T)
#map to human gene entrez id 
brain_sce = trans_mmu_to_hsa_stat(brain_sce, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
brain_sce_subclass = trans_mmu_to_hsa_stat(brain_sce_subclass, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
brain_sce_region_class = trans_mmu_to_hsa_stat(brain_sce_region_class, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
brain_sce_region_subclass = trans_mmu_to_hsa_stat(brain_sce_region_subclass, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
brain_sce_region_cluster = trans_mmu_to_hsa_stat(brain_sce_region_cluster, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")

#print magma files
print_magma_fuma_tbl(brain_sce, "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.fine_cluster.top10_magma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.fine_cluster.magma.aux.txt"))
print_magma_fuma_tbl(brain_sce_subclass, "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.subclass.top10_magma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.subclass.magma.aux.txt"))
print_magma_fuma_tbl(brain_sce_region_class, "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.region_class.top10_magma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.region_class.magma.aux.txt"))
print_magma_fuma_tbl(brain_sce_region_subclass, "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.region_subclass.top10_magma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.region_subclass.magma.aux.txt"))
print_magma_fuma_tbl(brain_sce_region_cluster, "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.region_cluster.top10_magma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.region_cluster.magma.aux.txt"))

###### fuma data preparation ##### 
#logcpm
brain_sce = scuttle::logNormCounts(brain_sce, assay.type = "cpm",size_factors=rep(1, ncol(brain_sce))) #no size factors
brain_sce_subclass = scuttle::logNormCounts(brain_sce_subclass, assay.type = "cpm",size_factors=rep(1, ncol(brain_sce_subclass))) 
brain_sce_region_class = scuttle::logNormCounts(brain_sce_region_class, assay.type = "cpm",size_factors=rep(1, ncol(brain_sce_region_class))) 
brain_sce_region_subclass = scuttle::logNormCounts(brain_sce_region_subclass, assay.type = "cpm",size_factors=rep(1, ncol(brain_sce_region_subclass))) 
brain_sce_region_cluster = scuttle::logNormCounts(brain_sce_region_cluster, assay.type = "cpm",size_factors=rep(1, ncol(brain_sce_region_cluster))) 

#calculate mean expression and map to human gene id (change assay)
brain_sce = cal_stat(brain_sce, meta_data = as.data.frame(colData(brain_sce)), group = "fine_cluster",assay_name ="logcounts" ,mean_only=T)
brain_sce_subclass = cal_stat(brain_sce_subclass, meta_data = as.data.frame(colData(brain_sce_subclass)), group = "subclass",assay_name ="logcounts" ,mean_only=T)
brain_sce_region_class = cal_stat(brain_sce_region_class, meta_data = as.data.frame(colData(brain_sce_region_class)), group = "region_class",assay_name ="logcounts" ,mean_only=T)
brain_sce_region_subclass = cal_stat(brain_sce_region_subclass, meta_data = as.data.frame(colData(brain_sce_region_subclass)), group = "region_subclass",assay_name ="logcounts" ,mean_only=T)
brain_sce_region_cluster = cal_stat(brain_sce_region_cluster, meta_data = as.data.frame(colData(brain_sce_region_cluster)), group = "region_cluster",assay_name ="logcounts" ,mean_only=T)
#map to entrez id
brain_sce = trans_mmu_to_hsa_stat(brain_sce, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
brain_sce_subclass = trans_mmu_to_hsa_stat(brain_sce_subclass, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
brain_sce_region_class = trans_mmu_to_hsa_stat(brain_sce_region_class, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
brain_sce_region_subclass = trans_mmu_to_hsa_stat(brain_sce_region_subclass, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")
brain_sce_region_cluster = trans_mmu_to_hsa_stat(brain_sce_region_cluster, gene_mapping_table = mmu_hsa_mapping, from = "mmu_symbol", to = "hsa_entrez")

#print fuma files
print_magma_fuma_tbl(brain_sce, "FUMA", main_table_path = here("data","expr","Saunders","Saunders.fine_cluster.fuma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.fine_cluster.fuma.aux.txt"))
print_magma_fuma_tbl(brain_sce_subclass, "FUMA", main_table_path = here("data","expr","Saunders","Saunders.subclass.fuma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.subclass.fuma.aux.txt"))
print_magma_fuma_tbl(brain_sce_region_class, "FUMA", main_table_path = here("data","expr","Saunders","Saunders.region_class.fuma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.region_class.fuma.aux.txt"))
print_magma_fuma_tbl(brain_sce_region_subclass, "FUMA", main_table_path = here("data","expr","Saunders","Saunders.region_subclass.fuma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.region_subclass.fuma.aux.txt"))
print_magma_fuma_tbl(brain_sce_region_cluster, "FUMA", main_table_path = here("data","expr","Saunders","Saunders.region_cluster.fuma.txt"),aux_table_path = here("data","expr","Saunders","Saunders.region_cluster.fuma.aux.txt"))

##scdrs 
#cells to keep (only cell types in the same analysis)
brain_sce_clean = brain_sce[,brain_sce$fine_cluster %in% names(get_seismic_ct_info(brain_sce, "cell_num"))]
brain_sce_subclass_clean = brain_sce_subclass[,brain_sce_subclass$subclass %in% names(get_seismic_ct_info(brain_sce_subclass, "cell_num"))]
brain_sce_region_class_clean = brain_sce_region_class[,brain_sce_region_class$region_class %in% names(get_seismic_ct_info(brain_sce_region_class, "cell_num"))]
brain_sce_region_subclass_clean = brain_sce_region_subclass[,brain_sce_region_subclass$region_subclass %in% names(get_seismic_ct_info(brain_sce_region_subclass, "cell_num"))]
brain_sce_region_cluster_clean = brain_sce_region_cluster[,brain_sce_region_cluster$region_cluster %in% names(get_seismic_ct_info(brain_sce_region_cluster, "cell_num"))]

#keep only anndata
brain_ad = AnnData(X = assay(brain_sce_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(brain_sce_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(brain_sce_clean)) %>% set_rownames(.$symbol) )
brain_subclass_ad = AnnData(X = assay(brain_sce_subclass_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_subclass_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                   obs = colData(brain_sce_subclass_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                   var = data.frame(symbol= rownames(brain_sce_subclass_clean)) %>% set_rownames(.$symbol) )
brain_region_class_ad = AnnData(X = assay(brain_sce_region_class_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_region_class_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                            obs = colData(brain_sce_region_class_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                            var = data.frame(symbol= rownames(brain_sce_region_class_clean)) %>% set_rownames(.$symbol) )
brain_region_subclass_ad = AnnData(X = assay(brain_sce_region_subclass_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_region_subclass_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                                obs = colData(brain_sce_region_subclass_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                                var = data.frame(symbol= rownames(brain_sce_region_subclass_clean)) %>% set_rownames(.$symbol) )
brain_region_cluster_ad = AnnData(X = assay(brain_sce_region_cluster_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_region_cluster_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                                obs = colData(brain_sce_region_cluster_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                                var = data.frame(symbol= rownames(brain_sce_region_cluster_clean)) %>% set_rownames(.$symbol) )

#write out h5ad
brain_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.fine_cluster.clean.h5ad"))
brain_subclass_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.subclass.clean.h5ad"))
brain_region_class_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.region_class.clean.h5ad"))
brain_region_subclass_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.region_subclass.clean.h5ad"))
brain_region_cluster_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.region_cluster.clean.h5ad"))
