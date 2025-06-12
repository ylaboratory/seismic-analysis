#prepare files for scDRS/FUMA/MAGMA 

##### 1. load packages and data#######
###load packages
if (!require("here")) {
  install.packages("here")
  library("here")
}

if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}

if (!require("SingleCellExperiment")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("SingleCellExperiment")
  library("SingleCellExperiment")
} 

#create anndata object and print out
if (!require("anndata")) {
  install.packages("anndata")
  library("anndata")
}


#load function
source(here("src","tools","magma_fuma_file_prep.R"))
source(here("src","tools","sparse_mat_util.R"))

#load the object that's filtered already in the Tabula_muris_analysis scripts
load(here("data","expr","Tabula_muris","TM_processed.rda"))
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
#calculate mean expression 
facs_mean <- calc_ct_mean(facs_obj_sce, assay_name = "counts", ct_label_col = "cluster_name")
droplet_mean <- calc_ct_mean(droplet_obj_sce, assay_name = "counts", ct_label_col = "cluster_name")

#filter genes without unique mapping and convert gene id
facs_mean <- facs_mean[, colnames(facs_mean) %in%mmu_hsa_mapping$mmu_symbol] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$mmu_symbol)])

droplet_mean <- droplet_mean[, colnames(droplet_mean) %in%mmu_hsa_mapping$mmu_symbol] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$mmu_symbol)])

#filter genes without expression
facs_mean <- facs_mean[, which(colSums(facs_mean)>0)]
droplet_mean = droplet_mean[, which(colSums(droplet_mean)>0)]

#convert to CPM
facs_mean <- sweep(facs_mean*1e6, MARGIN=1, STATS=rowSums(facs_mean), FUN="/")
droplet_mean <- sweep(droplet_mean*1e6, MARGIN=1, STATS=rowSums(droplet_mean), FUN="/")

#print magma files
print_magma_fuma_tbl(facs_mean, "MAGMA", main_table_path = here("data","expr","Tabula_muris","tm_facs.top10_magma.txt"),
                     aux_table_path = here("data","expr","Tabula_muris","tm_facs.magma.aux.txt"))
print_magma_fuma_tbl(droplet_mean, "MAGMA", main_table_path = here("data","expr","Tabula_muris","tm_droplet.top10_magma.txt"),
                     aux_table_path = here("data","expr","Tabula_muris","tm_droplet.magma.aux.txt"))


###### fuma data preparation ##### 
#cpm
assay(facs_obj_sce , "cpm") <- scuttle::calculateCPM(facs_obj_sce , assay.type = "counts", size.factors = colSums(assay(facs_obj_sce , "counts")))
assay(droplet_obj_sce, "cpm") <- scuttle::calculateCPM(droplet_obj_sce, assay.type = "counts", size.factors = colSums(assay(droplet_obj_sce, "counts")))
#logcpm
facs_obj_sce <- scuttle::logNormCounts(facs_obj_sce, assay.type = "cpm", size.factors=rep(1, ncol(facs_obj_sce)))  
droplet_obj_sce <- scuttle::logNormCounts(droplet_obj_sce, assay.type = "cpm",size.factors=rep(1,ncol(droplet_obj_sce))) 

#calculate mean expression and map to human gene id (change assay)
facs_mean_log <- calc_ct_mean(facs_obj_sce, assay_name = "logcounts", ct_label_col = "cluster_name")
droplet_mean_log <- calc_ct_mean(droplet_obj_sce, assay_name = "logcounts", ct_label_col = "cluster_name")

facs_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(facs_mean_log), from = "mmu_symbol")
droplet_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(droplet_mean_log), from = "mmu_symbol")

#print fuma files
print_magma_fuma_tbl(t(facs_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Tabula_muris","tm_facs.fuma.txt"),
                     aux_table_path = here("data","expr","Tabula_muris","tm_facs.fuma.aux.txt"))
print_magma_fuma_tbl(t(droplet_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Tabula_muris","tm_droplet.fuma.txt"),
                     aux_table_path = here("data","expr","Tabula_muris","tm_droplet.fuma.aux.txt"))

##scdrs 
#cells to keep (only cell types in the same analysis)
facs_obj_sce_clean <- facs_obj_sce[,facs_obj_sce$cluster_name %in% rownames(facs_mean)]
droplet_obj_sce_clean <- droplet_obj_sce[,droplet_obj_sce$cluster_name %in% rownames(droplet_mean)]

facs_ad = AnnData(X = assay(facs_obj_sce_clean,"counts")%>% set_rownames(rownames(rowData(facs_obj_sce_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(facs_obj_sce_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(facs_obj_sce_clean)) %>% set_rownames(.$symbol) )
droplet_ad = AnnData(X = assay(droplet_obj_sce_clean,"counts")%>% set_rownames(rownames(rowData(droplet_obj_sce_clean))) %>% set_colnames(paste0("cell",1:ncol(.))) %>% t(), 
                  obs = colData(droplet_obj_sce_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(droplet_obj_sce_clean)) %>% set_rownames(.$symbol) )

facs_ad %>%  write_h5ad(here("data","expr","Tabula_muris","facs.clean.h5ad"))
droplet_ad %>%  write_h5ad(here("data","expr","Tabula_muris","droplet.clean.h5ad"))



