#null simulation: take input data set
#1: expression data sets
#2: (cell) seed file, with each rows indicate the index of cells to be sampled as the fake cell type
#3: gs file
#4: gwas raw file
#5: process file: print out process of the program
#6: number of thread
#7: output file
#8: output directory: for temporary FUMA and MAGMA file
args <- commandArgs(trailingOnly = TRUE)
options(warn = -1)

#load packages
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("future"))
suppressMessages(library("future.apply"))
suppressMessages(library("here"))
suppressMessages(library("seismicGWAS"))
suppressMessages(library("SingleCellExperiment"))

#load source code
source(here("src","tools","magma_fuma_file_prep.R"))
source(here("src","tools","sparse_mat_util.R"))

#load data
data_path = args[1]
load(args[1])
cell_seed <- read.table(args[2], header=T,sep="\t")
z_score_file <- args[3]
trait_zscore <- read.table(z_score_file, header = T)
magma_raw_file <- args[4]
process_file <- args[5]
num_cores <- as.numeric(args[6])
output_file <- args[7]
temp_dir <- args[8]

print("NUM_CORE=",num_cores)

#data_path = "data/expr/null_sim/expr_rda_rs/expr_ds_1.rda"
#load(data_path)
#cell_seed = read.table("data/expr/null_sim/seed_table/sample_cell_idx.1.txt", header=T,sep="\t")
#z_score_file = "data/gwas/null_sim/zscore/gs_1.genes.out"
#trait_zscore = read.table("data/gwas/null_sim/zscore/gs_1.genes.out", header = T)
#magma_raw_file = "data/gwas/null_sim/magma_raw/gs_1.genes.raw"
#process_file = "data/test.txt"
#num_cores = 10
#output_file = "data/test_res.txt"
#=temp_dir = "data/temp"
#load gwas zscore and gene id mapping


#return a seismic p value for a target cell type
seismic_p_value <- function(data_sce, gwas_zscore_df, group, target_cell_type){
  seismic_sscore  <- calc_specificity(sce = data_sce , ct_label_col = group)
  seismic_sscore_hsa <- translate_gene_ids(seismic_sscore, from = "mmu_symbol")
  p_value_df <- get_ct_trait_associations(sscore = seismic_sscore_hsa, magma = gwas_zscore_df)
  return(p_value_df$pvalue[p_value_df$cell_type==target_cell_type])
}

#return a s-magma p value for target cell type
magma_p_value <- function(data_sce, temp_file_header, magma_raw_path, group, target_cell_type) {
  #prepare s-magma file
  mean_mat <- calc_ct_mean(data_sce, assay_name = "cpm", ct_label_col = group)
  mean_mat_hsa <- translate_gene_ids(t(mean_mat), from = "mmu_symbol")
  print_magma_fuma_tbl(t(mean_mat_hsa), "MAGMA", 
                       main_table_path = paste0(temp_file_header, ".magma.txt"),
                       aux_table_path = paste0(temp_file_header, ".magma.aux.txt"), verbose = F)
  
  #run magma
  ms <- system(paste0(here("bin", "magma", "magma"), " --gene-results ", magma_raw_path, " --set-annot ", temp_file_header, ".magma.txt --out ", temp_file_header,".magma"), intern = TRUE)
  
  #read in results
  magma_res_tbl <- read.table(paste0(temp_file_header, ".magma.gsa.out"), header = TRUE)
  magma_annot_tbl <- read.table(paste0(temp_file_header, ".magma.aux.txt"), header = TRUE, sep="\t")
  
  #get the encoded name of the cell type
  encoded_cell_type <- magma_annot_tbl$encoded_name[which(magma_annot_tbl$cell_type == target_cell_type)] 
  
  #get results
  res_p_value <- magma_res_tbl$P[magma_res_tbl$VARIABLE==encoded_cell_type]
  
  #remove unwanted files
  ms <- system(paste0("rm ", temp_file_header, ".*"), intern = TRUE)
  
  return(res_p_value)
}

fuma_p_value <- function(data_sce,  temp_file_header, magma_raw_path, group, target_cell_type) {
  #prepare fuma file
  mean_mat <- calc_ct_mean(data_sce, assay_name = "cpm", ct_label_col = group)
  mean_mat_hsa <- translate_gene_ids(t(mean_mat), from = "mmu_symbol")
  print_magma_fuma_tbl(t(mean_mat_hsa), "FUMA", 
                       main_table_path = paste0(temp_file_header, ".fuma.txt"),
                       aux_table_path = paste0(temp_file_header, ".fuma.aux.txt"), verbose = F)

  ms = system(paste0(here("bin", "magma", "magma"), " --gene-results ", magma_raw_path, " --gene-covar ", temp_file_header, ".fuma.txt --model condition-hide=Average direction=greater --out ", temp_file_header,".fuma"), intern = TRUE)
  
  #read in results
  fuma_res_tbl <- read.table(paste0(temp_file_header, ".fuma.gsa.out"), header = TRUE)
  fuma_annot_tbl <- read.table(paste0(temp_file_header, ".fuma.aux.txt"), header = TRUE, sep="\t")
  
  #get the encoded name of the cell type
  encoded_cell_type <- fuma_annot_tbl$encoded_name[which(fuma_annot_tbl$cell_type == target_cell_type)] 
  
  #get results
  res_p_value <- fuma_res_tbl$P[fuma_res_tbl$VARIABLE==encoded_cell_type]
  
  #remove unwanted files
  ms <- system(paste0("rm ", temp_file_header, ".*"), intern = TRUE)
  
  return(res_p_value)
}

#return a vector of p value
get_p_value <- function(data_sce, gwas_zscore_df, magma_raw_path, gene_mapping_table, group, cell_seed_vec, temp_file_header) {
  #set target cell type
  colData(data_sce)[[group]][cell_seed_vec] <- "fake_cell_type" 
  seismic_p <- seismic_p_value(data_sce, gwas_zscore_df = gwas_zscore_df, group = group, target_cell_type = "fake_cell_type")
  magma_p <- magma_p_value(data_sce,  temp_file_header = temp_file_header, magma_raw_path = magma_raw_path, group = group, target_cell_type = "fake_cell_type")
  fuma_p <- fuma_p_value(data_sce, temp_file_header = temp_file_header, magma_raw_path = magma_raw_path, group = group, target_cell_type = "fake_cell_type")
  
  #when FUMA return nothing because of correlation
  if (is.null(fuma_p)) {
    fuma_p = NA
  }
  
  #return
  return(c(seismic_p, magma_p, fuma_p))
}

#function wrapper
sim_p_value_by_idx <- function(i) {
  if (i < 1 | i > nrow(cell_seed)) {
    stop("i is not correct")
  }
  
  if (i %% 10 == 0) {
    cat(paste("Running seed ", i, "\n"), file = process_file, append = TRUE)
  }
  
  #get the target fake cell type index
  cell_seed_vec <- cell_seed[i, ] %>% unlist
  
  #get the rda name 
  data_name <- str_extract(data_path, pattern = "(?<=/)[^/]+$") %>% 
    gsub(pattern = ".rda", replacement = "", fixed = TRUE, x = .) 
  
  #p avlues
  p_val_vec = get_p_value(data_sce = sce, 
                          gwas_zscore_df = trait_zscore, 
                          group = "cell_ontology_class", 
                          cell_seed_vec = cell_seed_vec, 
                          temp_file_header = paste0(temp_dir, "/", data_name, ".", i),
                          magma_raw_path = magma_raw_file)
  
  return(p_val_vec)
}



library(parallel)
result = mclapply(1:10000, sim_p_value_by_idx, mc.cores = num_cores)

result_df = data.frame(index = 1:10000, seismic_p = unlist(map(result, ~.x[1])),magma_p = unlist(map(result, ~.x[2])),fuma_p = unlist(map(result, ~.x[3])))
write.table(result_df, file = output_file, sep = "\t",quote = F,row.names = F)
