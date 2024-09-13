#generate data
if (!require("here")){
  install.packages("here")
  library("here")
}
if (!require("magrittr")){
  install.packages("magrittr")
  library("magrittr")
}
if(!require("tidyverse")){
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("SingleCellExperiment")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SingleCellExperiment")
  library("SingleCellExperiment")
} 
library("scDesign3")
library("seismicGWAS")

#load objects
load(here("data","ref","mapping","mmu_hsa_mapping.rda"))

#load functions
source(here("src","tools","sparse_mat_util.R"))
source(here("src","tools","munge_sce_mat.R"))

#load("data/test.rda")
#load("data/test_new.rda")
#load("data/expr/causal_sim/ct_scdesign3_model/expr_ds_1.scdesign3_para.rda")
#load("data/expr/causal_sim/standard_scdesign3_model/expr_ds_1.scdesign3_para.rda")

##### function definition for scdDesign3 data generation ######
perturb_count_matrix <- function(sce, covar_df, cellids, geneids, effect_size_para, copula_para_list, para_list ){
  gene_idx <- match(geneids, rownames(sce))
  gene_loc <- match(geneids, colnames(copula_para_list$`1`)) #gene index in copula matrix
  cell_idx <- match(cellids, colnames(sce))
  covar_df <- dplyr::slice(covar_df, cell_idx)
  new_count_mat <- simu_new(
    sce = sce[gene_idx, cell_idx],
    mean_mat = para_list$mean_mat[cell_idx, gene_idx] %>% sweep(MARGIN = 1, STATS = effect_size_para, FUN = "*"),
    sigma_mat = para_list$sigma_mat[cell_idx, gene_idx],
    zero_mat = para_list$zero_mat[cell_idx, gene_idx],
    quantile_mat = NULL,
    copula_list = list("1" = copula_para_list$`1`[gene_loc, gene_loc]),
    n_cores = 10,
    family_use = "nb",
    input_data = mutate(covar_df, corr_group = 1),
    new_covariate = covar_df,
    important_feature = "all",
    filtered_gene = NULL #since all genes have expression in more than 2 cells
  )
  return(new_count_mat)
}



#### load data information
## load replace mat results
#parameter_df <- read.table("data/expr/causal_sim/sc3_gene_num_new/gene_400/perturbed_expr/parameter_df.txt", header = T) %>% as_tibble() %>% filter(expr_name == "expr_ds_1")
parameter_df <- read.table("data/expr/causal_sim/sc3_gene_num/gene_400/perturbed_expr/parameter_df.txt", header = T) %>% as_tibble() %>% filter(expr_name == "expr_ds_1")
i <- 1
ori_mat <- read.table(parameter_df$perturbed_mat_file[i], header = T) %>% as.matrix()
load(parameter_df$expr_rda_file[i])
data_sce <- sce
colnames(data_sce) <- data_sce$cellid

# add total counts
data_sce$library_size <- colSums(assay(data_sce, "counts"))

#extract cell information
gene_list <- rownames(ori_mat)
#cell_list <- colnames(ori_mat)
cell_list <- sample(colnames(data_sce), 100)
coldata_df <- colData(data_sce) %>% 
  as.data.frame() %>% 
  set_rownames(.$cellid) %>% 
  select(cell_ontology_class, library_size) 

#### load data parameter list by scDesign3 and generate new data
set.seed(456)
out_mat <- perturb_count_matrix(sce = data_sce, covar_df = coldata_df, cellids = cell_list, geneids = gene_list, 
                                effect_size_para = 2, copula_para_list = sce_copula$copula_list, para_list = sce_para)

#### process data
replace_mat <- out_mat
data_sce$cell_ontology_class[which(colnames(data_sce) %in% cell_list)] = "target_cell_type"

row_indices <- match(rownames(replace_mat), rownames(data_sce))
col_indices <- match(colnames(replace_mat), colnames(data_sce))

#modify count matrix
counts_mat <- as.matrix(assay(data_sce, "counts"))
size_factor <- colSums(counts_mat) / colSums(2^assay(data_sce, "logcounts") -1)  #size factor
counts_mat[row_indices, col_indices] <- as.matrix(replace_mat)
counts_mat <- as(counts_mat, "dgCMatrix")

logcounts_mat <- sweep_sparse(counts_mat, margin = 2, stats = size_factor) %>% transform_sparse(function(x) log2(x+1))
logcpm_mat <- sweep_sparse(counts_mat*1e6, margin = 2, stats = colSums(counts_mat)) %>% transform_sparse(function(x) log2(x+1))

assay(data_sce, "counts") <- counts_mat
assay(data_sce, "logcounts") <- logcounts_mat
assay(data_sce, "logcpm") <- logcpm_mat

#get association
sscore <- calc_specificity(data_sce, assay_name = "logcounts", ct_label_col = "cell_ontology_class")
sscore_hsa <- translate_gene_ids(sscore, from = "mmu_symbol")
all_association <- get_ct_trait_associations(sscore = sscore_hsa, magma = parameter_df$gs_zscore_file[i])

for (i in 1:10){
  ori_mat <- read.table(parameter_df$perturbed_mat_file[i], header = T) %>% as.matrix()
  load(parameter_df$expr_rda_file[i])
  data_sce <- sce
  colnames(data_sce) <- data_sce$cellid
  
  # add total counts
  data_sce$library_size <- colSums(assay(data_sce, "counts"))
  
  #extract cell information
  coldata_df <- colData(data_sce) %>% 
    as.data.frame() %>% 
    set_rownames(.$cellid) %>% 
    select(cell_ontology_class, library_size)
  
  gene_list <- rownames(ori_mat)
  cell_list <- colnames(ori_mat)
  
  #### load data parameter list by scDesign3 and generate new data
  set.seed(456)
  out_mat <- perturb_count_matrix(sce = data_sce, covar_df = coldata_df, cellids = cell_list, geneids = gene_list, 
                                  effect_size_para = 2, copula_para_list = sce_copula$copula_list, para_list = sce_para)
  
  #### process data
  replace_mat <- out_mat
  data_sce$cell_ontology_class[which(colnames(data_sce) %in% colnames(replace_mat))] = "target_cell_type"
  
  row_indices <- match(rownames(replace_mat), rownames(data_sce))
  col_indices <- match(colnames(replace_mat), colnames(data_sce))
  
  #modify count matrix
  counts_mat <- as.matrix(assay(data_sce, "counts"))
  size_factor <- colSums(counts_mat) / colSums(2^assay(data_sce, "logcounts") -1)  #size factor
  counts_mat[row_indices, col_indices] <- as.matrix(replace_mat)
  counts_mat <- as(counts_mat, "dgCMatrix")
  
  logcounts_mat <- sweep_sparse(counts_mat, margin = 2, stats = size_factor) %>% transform_sparse(function(x) log2(x+1))
  logcpm_mat <- sweep_sparse(counts_mat*1e6, margin = 2, stats = colSums(counts_mat)) %>% transform_sparse(function(x) log2(x+1))
  
  assay(data_sce, "counts") <- counts_mat
  assay(data_sce, "logcounts") <- logcounts_mat
  assay(data_sce, "logcpm") <- logcpm_mat
  
  #get association
  sscore <- calc_specificity(data_sce, assay_name = "logcounts", ct_label_col = "cell_ontology_class")
  sscore_hsa <- translate_gene_ids(sscore, from = "mmu_symbol")
  all_association <- get_ct_trait_associations(sscore = sscore_hsa, magma = parameter_df$gs_zscore_file[i])
  
  print(all_association$pvalue[all_association$cell_type == "target_cell_type"])
}
