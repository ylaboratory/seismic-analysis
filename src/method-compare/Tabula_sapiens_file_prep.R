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

if (!require("seismicGWAS")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

#create anndata object and print out
if (!require("anndata")) {
  install.packages("anndata")
  library("anndata")
}

#load function
source(here("src","tools","magma_fuma_file_prep.R"))
source(here("src","tools","sparse_mat_util.R"))

#load the object that's filtered already in the Tabula_sapeins_analysis scripts
load(here("data","expr","Tabula_sapiens","TS_processed.new.rda"))
load(here("data","ref","mapping","mmu_hsa_mapping.rda"))

#keep genes with only unique mapping
mmu_hsa_mapping <- mmu_hsa_mapping %>% 
  distinct(hsa_ensembl, hsa_entrez) %>%
  drop_na() %>%
  group_by(hsa_ensembl) %>% 
  filter(n()==1) %>%
  group_by(hsa_entrez) %>%
  filter(n()==1) %>%
  ungroup() %>%
  mutate(hsa_entrez = as.character(hsa_entrez))

###### magma data preparation #####
#calculate mean expression and map to human gene id
ts_mean <- calc_ct_mean(ts_obj, assay_name = "counts", ct_label_col = "cluster_name")

#filter genes without unique mapping and convert gene id
ts_mean = ts_mean[, colnames(ts_mean) %in%mmu_hsa_mapping$hsa_ensembl] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$hsa_ensembl)])

#filter genes without expression
ts_mean = ts_mean[, which(colSums(ts_mean)>0)]

#convert to CPM
ts_mean = sweep(ts_mean*1e6, MARGIN=1, STATS=rowSums(ts_mean), FUN="/")

#print magma files
print_magma_fuma_tbl(ts_mean, "MAGMA", main_table_path = here("data","expr","Tabula_sapiens","ts.top10_magma.new.txt"),
                     aux_table_path = here("data","expr","Tabula_sapiens","ts.magma.aux.new.txt"))


###### fuma data preparation ##### 
#cpm
assay(ts_obj, "cpm") <- scuttle::calculateCPM(ts_obj, assay.type = "counts", size.factors = colSums(assay(ts_obj, "counts")))

#logcpm
ts_obj <- scuttle::logNormCounts(ts_obj, assay.type = "cpm",size_factors=rep(1, ncol(ts_obj))) #no size factors

#calculate mean expression and map to human gene id (change assay)
ts_mean_log <- calc_ct_mean(ts_obj, assay_name = "logcounts", ct_label_col = "cluster_name")

ts_mean_log_hsa <- seismicGWAS::translate_gene_ids(t(ts_mean_log), from = "hsa_ensembl")

#print fuma files
print_magma_fuma_tbl(t(ts_mean_log_hsa), "FUMA", main_table_path = here("data","expr","Tabula_sapiens","ts.fuma.new.txt"),
                     aux_table_path = here("data","expr","Tabula_sapiens","ts.fuma.aux.new.txt"))


##scdrs 
#cells to keep (only cell types in the same analysis)
ts_obj_clean <- ts_obj[,ts_obj$cluster_name %in% rownames(ts_mean)]

ts_count_mat <- assay(ts_obj_clean,"counts")%>% 
  set_rownames(rownames(rowData(ts_obj_clean))) %>% 
  set_colnames(paste0("cell",1:ncol(.)))

ts_count_mat_symbol <- seismicGWAS::translate_gene_ids(ts_count_mat, from = "hsa_ensembl", to = "hsa_symbol", multi_mapping = "sum")

ts_ad <- AnnData(X = t(ts_count_mat_symbol), 
                  obs = colData(ts_obj_clean) %>% as.data.frame %>% set_rownames(paste0("cell",1:nrow(.))), 
                  var = data.frame(symbol= rownames(ts_count_mat_symbol)) %>% set_rownames(.$symbol) )

ts_ad %>%  write_h5ad(here("data","expr","Tabula_sapiens","ts.clean.new.h5ad"))

