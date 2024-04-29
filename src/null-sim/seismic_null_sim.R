#null simulation: take input data set
#1: expression data sets
#2: (cell) seed file, with each rows indicate the index of cells to be sampled as the fake cell type
#3: gs file
#4: gwas raw file
#5: process file: print out process of the program
#6: number of thread
#7: output file
#8: output directory: for temporary FUMA and MAGMA file
args = commandArgs(trailingOnly = TRUE)
options(warn = -1)

#load packages
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("future"))
suppressMessages(library("future.apply"))
suppressMessages(library("here"))
suppressMessages(library("seismicGWAS"))
suppressMessages(library("SingleCellExperiment"))

#load data
data_path = args[1]
load(args[1])
cell_seed = read.table(args[2], header=T,sep="\t")
z_score_file = args[3]
trait_zscore <- read.table(z_score_file, header = T)
process_file = args[4]
num_cores= as.numeric(args[5])
output_file = args[6]

#return a seismic p value for a target cell type
seismic_p_value <- function(data_sce, gwas_zscore_df, group, target_cell_type){
  seismic_sscore  <- calc_specificity(sce = data_sce , ct_label_col = group)
  seismic_sscore_hsa <- translate_gene_ids(seismic_sscore, from = "mmu_symbol")
  p_value_df <- get_ct_trait_associations(sscore = seismic_sscore_hsa, magma = gwas_zscore_df)
  return(p_value_df$pvalue[p_value_df$cell_type==target_cell_type])
}

get_p_value <- function(data_sce, gwas_zscore_df, group, cell_seed_vec, temp_file_header) {
  #set target cell type
  colData(data_sce)[[group]][cell_seed_vec] <- "fake_cell_type" 
  seismic_p <- seismic_p_value(data_sce, gwas_zscore_df = gwas_zscore_df, group = group, target_cell_type = "fake_cell_type")
  
  return(seismic_p)
}


#function wrapper
sim_p_value_by_idx = function(i){
  if( i < 1 | i > nrow(cell_seed)){
    stop("i is not correct")
  }
  if(i %% 10 == 0) {
    cat(paste("Running seed ", i, "\n"), file = process_file, append = TRUE)
  }
  #get the target fake cell type index
  cell_seed_vec <- cell_seed[i, ] %>% unlist
  
  #get the rda name 
  data_name <- str_extract(data_path, pattern = "(?<=/)[^/]+$") %>% 
    gsub(pattern = ".rda", replacement = "", fixed = TRUE, x = .) 
  
  p_val = get_p_value(data_sce = sce , gwas_zscore_df = trait_zscore, 
                       group = "cell_ontology_class",  cell_seed_vec = cell_seed_vec)
  return(p_val)
}

#plan(multisession, workers = num_cores)
#result = future_sapply(1:100,sim_p_value_by_idx)

library(parallel)
result = mclapply(1:10000, sim_p_value_by_idx, mc.cores = num_cores)

result_df = data.frame(index = 1:10000, ours_p = unlist(result))
write.table(result_df, file = output_file, sep = "\t",quote = F,row.names = F)
