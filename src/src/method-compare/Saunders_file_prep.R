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
load(here("data","ref","mapping","mmu_hsa_mapping.rda"))

#keep genes with only unique mapping
mmu_hsa_mapping <- mmu_hsa_mapping %>% 
  distinct(mmu_symbol, hsa_entrez) %>%
  drop_na() %>%
  group_by(mmu_symbol) %>% 
  filter(n()==1) %>%
  group_by(hsa_entrez) %>%
  filter(n()==1) %>%
  ungroup() %>%
  mutate(hsa_entrez = as.character(hsa_entrez))

###### magma data preparation #####
#calculate mean expression and map to human gene id
fine_cluster_mean <- calc_ct_mean(brain_sce, assay_name = "counts", ct_label_col = "fine_cluster")
subclass_mean <- calc_ct_mean(brain_sce, assay_name = "counts", ct_label_col = "subclass")
region_subclass_mean <- calc_ct_mean(brain_sce, assay_name = "counts", ct_label_col = "region_subclass")
region_class_mean <- calc_ct_mean(brain_sce, assay_name = "counts", ct_label_col = "region_class")
region_cluster_mean <- calc_ct_mean(brain_sce, assay_name = "counts", ct_label_col = "region_cluster")

#filter genes without unique mapping and convert gene id
fine_cluster_mean <- fine_cluster_mean[, colnames(fine_cluster_mean) %in%mmu_hsa_mapping$mmu_symbol] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$mmu_symbol)])

subclass_mean <- subclass_mean[, colnames(subclass_mean) %in%mmu_hsa_mapping$mmu_symbol] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$mmu_symbol)])

region_subclass_mean <- region_subclass_mean[, colnames(region_subclass_mean) %in%mmu_hsa_mapping$mmu_symbol] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$mmu_symbol)])

region_class_mean <- region_class_mean[, colnames(region_class_mean) %in%mmu_hsa_mapping$mmu_symbol] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$mmu_symbol)])

region_cluster_mean <- region_cluster_mean[, colnames(region_cluster_mean) %in%mmu_hsa_mapping$mmu_symbol] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$mmu_symbol)])

#filter genes without expression
fine_cluster_mean <- fine_cluster_mean[, which(colSums(fine_cluster_mean)>0)]
subclass_mean <- subclass_mean[, which(colSums(subclass_mean)>0)]
region_subclass_mean <- region_subclass_mean[, which(colSums(region_subclass_mean)>0)]
region_class_mean <- region_class_mean[, which(colSums(region_class_mean)>0)]
region_cluster_mean <- region_cluster_mean[, which(colSums(region_cluster_mean)>0)]

#convert to CPM
fine_cluster_mean <- sweep(fine_cluster_mean*1e6, MARGIN=1, STATS=rowSums(fine_cluster_mean), FUN="/")
subclass_mean <- sweep(subclass_mean*1e6, MARGIN=1, STATS=rowSums(subclass_mean), FUN="/")
region_subclass_mean <- sweep(region_subclass_mean*1e6, MARGIN=1, STATS=rowSums(region_subclass_mean), FUN="/")
region_class_mean <- sweep(region_class_mean*1e6, MARGIN=1, STATS=rowSums(region_class_mean), FUN="/")
region_cluster_mean <- sweep(region_cluster_mean*1e6, MARGIN=1, STATS=rowSums(region_cluster_mean), FUN="/")

#print magma files
print_magma_fuma_tbl(fine_cluster_mean , "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.fine_cluster.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.fine_cluster.magma.aux.txt"))
print_magma_fuma_tbl(subclass_mean, "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.subclass.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.subclass.magma.aux.txt"))
print_magma_fuma_tbl(region_subclass_mean, "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.region_subclass.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.region_subclass.magma.aux.txt"))
print_magma_fuma_tbl(region_class_mean, "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.region_class.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.region_class.magma.aux.txt"))
print_magma_fuma_tbl(region_cluster_mean , "MAGMA", main_table_path = here("data","expr","Saunders","Saunders.region_cluster.top10_magma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.region_cluster.magma.aux.txt"))

###### fuma data preparation ##### 
#cpm
assay(brain_sce, "cpm")  <- scuttle::calculateCPM(brain_sce, assay.type = "counts", size.factors = colSums(assay(brain_sce, "counts")))

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
print_magma_fuma_tbl(t(fine_cluster_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","Saunders.fine_cluster.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.fine_cluster.fuma.aux.txt"))
print_magma_fuma_tbl(t(subclass_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","Saunders.subclass.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.subclass.fuma.aux.txt"))
print_magma_fuma_tbl(t(region_subclass_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","Saunders.region_subclass.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.region_subclass.fuma.aux.txt"))
print_magma_fuma_tbl(t(region_class_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","Saunders.region_class.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.region_class.fuma.aux.txt"))
print_magma_fuma_tbl(t(region_cluster_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Saunders","Saunders.region_cluster.fuma.txt"),
                     aux_table_path = here("data","expr","Saunders","Saunders.region_cluster.fuma.aux.txt"))

##scdrs 
#cells to keep (only cell types in the same analysis)
brain_sce_clean <- brain_sce[,brain_sce$fine_cluster %in% rownames(fine_cluster_mean)]
brain_sce_subclass_clean <- brain_sce[,brain_sce$subclass %in% rownames(subclass_mean)]
brain_sce_region_class_clean <- brain_sce[,brain_sce$region_class %in% rowames(region_class_mean)]
brain_sce_region_subclass_clean <- brain_sce[,brain_sce$region_subclass %in% rownames(region_subclass_mean)]
brain_sce_region_cluster_clean <- brain_sce[,brain_sce$region_cluster %in% rownames(region_cluster_mean)]

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
