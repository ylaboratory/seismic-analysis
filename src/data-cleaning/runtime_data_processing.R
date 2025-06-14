if (!require("reshape2")) {
  install.packages("reshape2")
  library("reshape2")
}

if (!require("here")) {
  install.packages("here")
  library("here")
}

if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}

if (!require("scran")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("scran")
  library("scran")
}

if (!require("zellkonverter")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("zellkonverter")
  library("zellkonverter")
} 
if(!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}
library(anndata)

#load function
source(here("src","tools","munge_sce_mat.R"))

# load data and filter data
ts_obj <- readH5AD(here("raw","expr","Tabula_sapiens","TabulaSapiens.h5ad"),reader = "R")

# quality control and filter cells 
ts_obj$mt_counts <- assay(ts_obj, "decontXcounts")[grep(pattern = "^MT-", rowData(ts_obj)$gene_symbol), ] %>% colSums()
ts_obj$mt_ratio <- ts_obj$mt_counts/ts_obj$n_counts_UMIs 

ts_obj <- ts_obj[, ts_obj$mt_ratio<=0.15] #keep more cells
ts_obj <- ts_obj[, ts_obj$method=="10X"] #most of the cells

# filter genes
ts_obj <- ts_obj[( rowSums( assay(ts_obj, "decontXcounts"))>=20),]
ts_obj <- ts_obj %>%
  set_rownames(rowData(ts_obj)$ensemblid %>% strsplit(split=".",fixed=T) %>% map(~.[1]) %>% unlist) %>%
  set_colnames(colData(ts_obj)$cell_id)

# normalization
cluster <- quickCluster(assay(ts_obj, "decontXcounts")) 
size.factor <- calculateSumFactors(assay(ts_obj, "decontXcounts"), cluster=cluster, min.mean=0.1)
ts_obj <- logNormCounts(ts_obj, size.factors = size.factor, assay.type = "decontXcounts" )
assays(ts_obj) <- assays(ts_obj)[c("decontXcounts","logcounts")]
ts_obj_symbol <- munge_sce_mat(ts_obj, rowData(ts_obj) %>% as_tibble(rownames = "ensembl") %>% select(ensembl, gene_symbol))

# subset data set
all_k <- rev(c(300000, 250000, 200000, 150000,100000, 50000,25000, 10000))
set.seed(101)
for (i in 1:5) {
  for (k in all_k) {
    seed_sample <- sample(ncol(ts_obj), k)
    seed_table <- data.frame(cell = paste0("cell_",1:k), cell_idx = seed_sample)
    write.table(seed_table, file=here("data","expr","seed_table",paste0("ds_",i), paste0("seed_table.",k/1000,"k.txt")), quote = F, col.names = T, row.names = T, sep="\t")
    #sce for ours
    ts.sce <- ts_obj[,seed_sample]
    save(ts.sce, file = here("data","expr","runtime","expr_rda",paste0("ds_",i), paste0("sample.",k/1000,"k.rda") ) )
    #sce for anndata
    ts.adata <- ts_obj_symbol[,seed_sample]
    ts.adata <- SCE2AnnData(ts.adata)
    ts.adata %>% write_h5ad(here("data","expr","runtime","expr_h5ad",paste0("ds_",i), paste0("sample.",k/1000,"k.h5ad")))
  }
}


# upsampling 
extra_k <- c(400000, 500000)
set.seed(100)
for (i in 1:5) {
  for (k in extra_k) {
    #sample seed
    seed_sample <- sample(ncol(ts_obj), k, replace = T)
    seed_table <- data.frame(cell = paste0("cell_",1:k), cell_idx = seed_sample)
    write.table(seed_table, file=here("data","expr","seed_table",paste0("ds_",i), paste0("seed_table.",k/1000,"k.txt")), quote = F, col.names = T, row.names = T, sep="\t")
    
    #sce for ours
    ts.sce <- ts_obj[,seed_sample]
    
    colData(ts.sce) <- colData(ts.sce) %>% 
      as.data.frame() %>%
      select(cell_id, cell_ontology_class, n_counts_UMIs, sizeFactor) %>%
      group_by(cell_id) %>%
      mutate(order = 1:n()) %>%
      ungroup() %>%
      mutate(cell_id = ifelse(order > 1, paste0(cell_id,":",order), cell_id)) %>%
      mutate(cell_ontology_class = ifelse(order > 1, paste0(cell_ontology_class, ":", order), cell_ontology_class)) %>%
      DataFrame() %>%
      set_rownames(.$cell_id)
    
    colnames(ts.sce) <- ts.sce$cell_id
    
    save(ts.sce, file = here("data","expr","runtime","expr_rda",paste0("ds_",i), paste0("sample.",k/1000,"k.rda")))
    
    #sce for anndata
    ts.adata <- ts_obj_symbol[,seed_sample]
    
    colData(ts.adata) <- colData(ts.sce)
    
    colnames(ts.adata) <- ts.adata$cell_id
    
    ts.adata <- SCE2AnnData(ts.adata)
    
    ts.adata %>% write_h5ad(here("data","expr","runtime","expr_h5ad",paste0("ds_",i), paste0("sample.",k/1000,"k.h5ad"))) #write out data
  }
}

