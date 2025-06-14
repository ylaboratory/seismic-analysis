# script for null simulation
if (!require("here")) {
  install.packages("here")
  library("here")
}
if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}
if(!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("scran")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("scran")
  library("scran")
}
if (!require("anndata")) {
  install.packages("anndata")
  library("anndata")
}

load(here("data","ref","mapping","mmu_hsa_mapping.rda"))
source(here("src","tools","sparse_mat_util.R"))

# sample gene sets and generate random trait Z-score
z_score_file <- list.files(here("data","gwas","tm_gwas","zscore"), full.names = T) %>% 
  set_names(str_extract(. , pattern = "(?<=/)[^/]+$" ) %>% gsub(pattern = ".genes.out",replacement = "", fixed = T )) %>%
  .[match(unique(names(.)), names(.))]

magma_raw_file <- list.files(here("data","gwas","tm_gwas","magma_raw"), full.names = T) %>% 
  set_names(str_extract(. , pattern = "(?<=/)[^/]+$" ) %>% gsub(pattern = ".genes.raw",replacement = "", fixed = T )) %>%
  .[match(unique(names(.)), names(.))]

z_score <- z_score_file %>% 
  map(~read.table(.x, header=T) %>% as_tibble()) %>%
  map(~select(.x,contains(c("GENE","ZSTAT")))) 


# sample gene sets
set.seed(98)
selected_traits <- sample(names(z_score_file),10)

# write out the sampled gene sets
map2(selected_traits, 1:10 , ~system(paste0("cp ",z_score_file[[.x]]," ",here("data","null_sim","zscore",paste0("gs_",.y,".genes.out")))))
map2(selected_traits, 1:10 , ~system(paste0("cp ",magma_raw_file[[.x]]," ",here("data","null_sim","magma_raw",paste0("gs_",.y,".genes.raw")))))

#write out the shuffled gene sets 
set.seed(100)
for (i in 1:10) {
  #selected trait
  trait = selected_traits[i]
  
  #select the corresonding data
  z_score_vec <- z_score[[trait]] %>% 
    set_colnames(c("GENE", "Trait")) 
  
  magma_raw_value <- magma_raw[[trait]]
  
  #print scdrs files
  z_score_vec %>% 
    left_join(mmu_hsa_mapping %>% distinct(hsa_entrez,hsa_symbol), by = c("GENE"="hsa_entrez"),multiple = "all") %>%
    relocate(hsa_symbol) %>%
    distinct() %>%
    drop_na(hsa_symbol) %>%
    select(-GENE) %>% 
    set_colnames(c("GENE",paste0("Trait_",i))) %>%
    group_by(GENE) %>% #get average by genes 
    mutate_if(is.numeric, ~mean(.)) %>%
    ungroup %>% 
    distinct() %>% 
    write_tsv(here("data","null_sim","scdrs_gs",paste0("gs_",i,".tsv")))
  
  #random shuffle 
  shuffle_idx = sample(1:nrow(z_score_vec))
  z_score_vec_rs <- z_score_vec %>% mutate(GENE = GENE[shuffle_idx]) #random shuffle gene index
  
  #print scdrs files
  z_score_vec_rs %>%
    left_join(mmu_hsa_mapping %>% distinct(hsa_entrez,hsa_symbol), by = c("GENE"="hsa_entrez"),multiple = "all") %>%
    relocate(hsa_symbol) %>%
    distinct() %>%
    drop_na(hsa_symbol) %>%
    select(-GENE) %>% 
    set_colnames(c("GENE",paste0("Trait_",i))) %>%
    group_by(GENE) %>% #get average by genes 
    mutate_if(is.numeric, ~mean(.)) %>%
    ungroup %>% 
    distinct() %>% 
    write_tsv(here("data","null_sim","scdrs_gs_rs",paste0("gs_",i,".tsv")))
  
  #write out the shuffled z-score vector
  z_score_vec_rs %>% 
    set_colnames(c("GENE","ZSTAT")) %>% 
    write.table(file= here("data","null_sim","zscore_rs",paste0("gs_",i,".genes.out")),col.names = T, row.names = F, quote = F )
}

# random expression data generation
load(here("data","expr","Tabula_muris","facs_clean.rda"))

#remove these witout cell ontology id: because it will be used as later analysis granularity
facs_obj <- facs_obj_sce[,which(!is.na(facs_obj_sce$cell_ontology_class))]

# normalization
cluster.facs <- quickCluster(facs_obj_sce, assay.type = "counts") 
facs.factor <- calculateSumFactors(facs_obj_sce, cluster=cluster.facs, min.mean=0.1)
facs_obj_sce <- logNormCounts(facs_obj_sce, size.factors = facs.factor )

# adjusted size factor: to make sure that logcounts and logCPM assays coexist
facs.adjusted.log <- assay(facs_obj_sce, "logcounts") %>% 
  as("dgCMatrix") 

# extract other matrix with different normalization 
facs.count_mat <-  assay(facs_obj_sce, "counts") %>%  
  as("dgCMatrix") #for scdrs

#cpm
facs.cpm <- facs.count_mat %>% 
  sweep_sparse(margin=2,stats = 1e6/colSums(facs.count_mat),fun="*") #for MAGMA

#logcpm
facs.logcpm <- facs.cpm %>% 
  sweep_sparse(margin = 2,stats = 1, fun = "+") %>%
  transform_sparse(fun = "log2") #for FUMA

# random sample and output 
set.seed(99)
for (i in 1:10) {
  # random select by cell index
  cell_idx <- sample(1:ncol(facs.adjusted.log), size=10000) #select 10000 cells
  
  # random select gene indx
  gene_idx <- sample(1:nrow(facs.adjusted.log)) #gene index
  
  # subset data set
  cell_anno <- colData(facs_obj_sce) %>% 
    as_tibble() %>% 
    mutate(cellid = paste0("cell_",1:n())) %>%
    dplyr::slice(cell_idx) %>%
    select(cellid, tissue, cell_ontology_class) %>%
    set_colnames(c("cellid","tissue","cell_ontology_class"))
  
  # write out for our pipeline
  sce <- SingleCellExperiment::SingleCellExperiment(list(logcounts = facs.adjusted.log[gene_idx,cell_idx] , 
                                                        cpm= facs.cpm[gene_idx, cell_idx],
                                                        logcpm = facs.logcpm[gene_idx,cell_idx]),
                                                   colData=DataFrame(cell_anno))
  
  #save expression dataset with random gene symbol permutation
  #i.e. reset the gene symbols to the original 
  sce <- sce %>% set_rownames(rownames(facs.adjusted.log)) 
  
  save(sce,cell_anno, file = here("data","null_sim","expr_rda_rs",paste0("expr_ds_",i,".rda")) )
  
  #write out anndata
  ann_out <- AnnData(X= t(facs.count_mat[gene_idx, cell_idx]),
                    obs=cell_anno %>% as.data.frame() %>% set_rownames(.$cellid),
                    var =  data.frame(symbol= rownames(facs.count_mat)[gene_idx]) %>% set_rownames(.$symbol))
  
  ann_out_rs <- AnnData(X = t(facs.count_mat[gene_idx, cell_idx]) %>% set_colnames(rownames(facs.count_mat)),
                       obs=cell_anno %>% as.data.frame() %>% set_rownames(.$cellid),
                       var = data.frame(symbol= rownames(facs.count_mat)) %>% set_rownames(.$symbol))
  
  ann_out %>% write_h5ad(here("data","null_sim","expr_h5ad",paste0("expr_ds_",i,".h5ad")))
  
  ann_out_rs %>% write_h5ad(here("data","null_sim","expr_h5ad_rs",paste0("expr_ds_",i,".h5ad")))
}


# random sample seeds (cell index) for later analysis: simu_data_gene.R
set.seed(101)
for (i in 1:10) {
  idx_mat <- matrix(nrow=10000,ncol=100)
  #index 
  for (j in 1:10000) {
    idx_mat[j,] <- sample(1:10000,size=100)
  }
  idx_mat <- idx_mat %>% 
    set_rownames(paste0("draw_",1:10000)) %>%
    set_colnames(paste0("cell_",1:100))
  idx_mat %>%  write.table( here("data","expr","null_sim","seed_table",paste0("sample_cell_idx.",i,".txt")),sep="\t",col.names = T, row.names = T, quote = F)
}
