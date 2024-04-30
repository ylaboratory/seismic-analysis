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
load(here("data","expr","Saunders","Saunders_processed.rda"))

###### magma data preparation #####
#cpm
assay(brain_sce, "cpm")  <- scuttle::calculateCPM(brain_sce, assay.type = "counts", size.factors = colSums(assay(brain_sce, "counts")))

#calculate mean expression and map to human gene id
fine_cluster_mean <- calc_ct_mean(brain_sce, assay_name = "cpm", ct_label_col = "fine_cluster")
subclass_mean <- calc_ct_mean(brain_sce, assay_name = "cpm", ct_label_col = "subclass")
region_subclass_mean <- calc_ct_mean(brain_sce, assay_name = "cpm", ct_label_col = "region_subclass")
region_class_mean <- calc_ct_mean(brain_sce, assay_name = "cpm", ct_label_col = "region_class")
region_cluster_mean <- calc_ct_mean(brain_sce, assay_name = "cpm", ct_label_col = "region_cluster")

#map to human gene entrez id 
fine_cluster_mean_hsa <- seismicGWAS::translate_gene_ids(t(fine_cluster_mean), from = "mmu_symbol")
subclass_mean_hsa <- seismicGWAS::translate_gene_ids(t(subclass_mean), from = "mmu_symbol")
region_subclass_mean_hsa <- seismicGWAS::translate_gene_ids(t(region_subclass_mean), from = "mmu_symbol")
region_class_mean_hsa <- seismicGWAS::translate_gene_ids(t(region_class_mean), from = "mmu_symbol")
region_cluster_mean_hsa <- seismicGWAS::translate_gene_ids(t(region_cluster_mean), from = "mmu_symbol")

#print magma files
print_magma_fuma_tbl(t(fine_cluster_mean_hsa) , "MAGMA", main_table_path = here("data","expr","Saunders","new_Saunders.fine_cluster.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.fine_cluster.magma.aux.txt"))
print_magma_fuma_tbl(t(subclass_mean_hsa), "MAGMA", main_table_path = here("data","expr","Saunders","new_Saunders.subclass.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.subclass.magma.aux.txt"))
print_magma_fuma_tbl(t(region_subclass_mean_hsa), "MAGMA", main_table_path = here("data","expr","Saunders","new_Saunders.region_class.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.region_class.magma.aux.txt"))
print_magma_fuma_tbl(t(region_class_mean_hsa), "MAGMA", main_table_path = here("data","expr","Saunders","new_Saunders.region_subclass.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.region_subclass.magma.aux.txt"))
print_magma_fuma_tbl(t(region_cluster_mean_hsa) , "MAGMA", main_table_path = here("data","expr","Saunders","new_Saunders.region_cluster.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.region_cluster.magma.aux.txt"))

###### fuma data preparation ##### 
#logcpm
brain_sce <- scuttle::logNormCounts(brain_sce, assay.type = "cpm",size_factors=rep(1, ncol(brain_sce))) #no size factors

#calculate mean expression and map to human gene id
fine_cluster_mean_log <- calc_ct_mean(brain_sce, assay_name = "logcounts", ct_label_col = "fine_cluster")
subclass_mean_log <- calc_ct_mean(brain_sce, assay_name = "logcounts", ct_label_col = "subclass")
region_subclass_mean_log <- calc_ct_mean(brain_sce, assay_name = "logcounts", ct_label_col = "region_subclass")
region_class_mean_log <- calc_ct_mean(brain_sce, assay_name = "logcounts", ct_label_col = "region_class")
region_cluster_mean_log <- calc_ct_mean(brain_sce, assay_name = "logcounts", ct_label_col = "region_cluster")

#map to human gene entrez id 
fine_cluster_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(fine_cluster_mean_log), from = "mmu_symbol")
subclass_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(subclass_mean_log), from = "mmu_symbol")
region_subclass_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(region_subclass_mean_log), from = "mmu_symbol")
region_class_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(region_class_mean_log), from = "mmu_symbol")
region_cluster_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(region_cluster_mean_log), from = "mmu_symbol")

#print fuma files
print_magma_fuma_tbl(t(fine_cluster_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","new_Saunders.fine_cluster.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.fine_cluster.fuma.aux.txt"))
print_magma_fuma_tbl(t(subclass_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","new_Saunders.subclass.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.subclass.fuma.aux.txt"))
print_magma_fuma_tbl(t(region_subclass_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","new_Saunders.region_subclass.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.region_subclass.fuma.aux.txt"))
print_magma_fuma_tbl(t(region_class_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","new_Saunders.region_class.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.region_class.fuma.aux.txt"))
print_magma_fuma_tbl(t(region_cluster_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","new_Saunders.region_cluster.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","new_Saunders.region_cluster.fuma.aux.txt"))

##scdrs 
#cells to keep (only cell types in the same analysis)
brain_sce_clean <- brain_sce[,brain_sce$fine_cluster %in% colnames(fine_cluster_mean_hsa)]
brain_sce_subclass_clean <- brain_sce[,brain_sce$subclass %in% colnames(subclass_mean_hsa)]
brain_sce_region_class_clean <- brain_sce[,brain_sce$region_class %in% colnames(region_class_mean_hsa)]
brain_sce_region_subclass_clean <- brain_sce[,brain_sce$region_subclass %in% colnames(region_subclass_mean_hsa)]
brain_sce_region_cluster_clean <- brain_sce[,brain_sce$region_cluster %in% colnames(region_cluster_mean_hsa)]

#keep only anndata
brain_ad <- AnnData(X = assay(brain_sce_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(brain_sce_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(brain_sce_clean)) %>% set_rownames(.$symbol) )
brain_subclass_ad <- AnnData(X = assay(brain_sce_subclass_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_subclass_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                   obs = colData(brain_sce_subclass_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                   var = data.frame(symbol= rownames(brain_sce_subclass_clean)) %>% set_rownames(.$symbol) )
brain_region_class_ad <- AnnData(X = assay(brain_sce_region_class_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_region_class_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                            obs = colData(brain_sce_region_class_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                            var = data.frame(symbol= rownames(brain_sce_region_class_clean)) %>% set_rownames(.$symbol) )
brain_region_subclass_ad <- AnnData(X = assay(brain_sce_region_subclass_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_region_subclass_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                                obs = colData(brain_sce_region_subclass_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                                var = data.frame(symbol= rownames(brain_sce_region_subclass_clean)) %>% set_rownames(.$symbol) )
brain_region_cluster_ad <- AnnData(X = assay(brain_sce_region_cluster_clean,"counts")%>% set_rownames(rownames(rowData(brain_sce_region_cluster_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                                obs = colData(brain_sce_region_cluster_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                                var = data.frame(symbol= rownames(brain_sce_region_cluster_clean)) %>% set_rownames(.$symbol) )

#write out h5ad
brain_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.fine_cluster.clean.h5ad"))
brain_subclass_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.subclass.clean.h5ad"))
brain_region_class_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.region_class.clean.h5ad"))
brain_region_subclass_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.region_subclass.clean.h5ad"))
brain_region_cluster_ad %>%  write_h5ad(here("data","expr","Saunders","Saunders.region_cluster.clean.h5ad"))
