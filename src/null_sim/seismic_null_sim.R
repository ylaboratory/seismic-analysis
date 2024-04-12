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
load(args[1])
cell_seed = read.table(args[2], header=T,sep="\t")
z_score_file = args[3]
#magma_raw_file = args[4]
process_file = args[4]
num_cores= as.numeric(args[5])
output_file = args[6]
#temp_dir = args[8]

#load gwas zscore and gene id mapping
trait_zscore = load_zscore(z_score_file)
data("mmu_hsa_mapping")  

#return a data_sce
seismic_gwas_p_value = function(data_sce, gwas_zscore, gene_mapping_table, group){ #retreive association of #fake_type
  data_sce = cal_stat(data_obj = data_sce, meta_data = as.data.frame(colData(data_sce)), group = group, assay_name="logcounts") #this is log normalized counts 
  data_sce = cal_sscore(data_obj = data_sce) 
  data_sce = trans_mmu_to_hsa_stat(data_sce, gene_mapping_table=gene_mapping_table, from="mmu_symbol", to="hsa_entrez")
  data_sce = add_glob_stats(data_sce, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 
  data_sce = cal_ct_asso(data_sce, gwas_zscore, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.1& max_exp_ct>0.1")
  p_value = get_ct_asso(data_sce, trait_name = "all",asso_model = "linear")[[1]]  %>% filter(cell_type == "fake_cell_type") %>% pull(Pvalue)
  return(p_value)
}


get_p_value = function(data_sce, gwas_zscore, gene_mapping_table, group, cell_seed_vec){
  metadata(data_sce) = list(NULL)
  colData(data_sce)[[group]][cell_seed_vec] = "fake_cell_type" #fake cell type
  our_p = seismic_gwas_p_value(data_sce, gwas_zscore = gwas_zscore, gene_mapping_table = gene_mapping_table, group = group)
  #return(c(our_p,magma_p,fuma_p))
  return(our_p)
}


#function wrapper
sim_p_value_by_idx = function(i){
  if( i < 1 | i > nrow(cell_seed)){
    stop("i is not correct")
  }
  if(i %% 10 == 0) {
    cat(paste("Running seed ", i, "\n"), file = process_file, append = TRUE)
  }
  cell_seed_vec = cell_seed[i,] %>% unlist
  #data_name = str_extract(args[1],pattern = "(?<=/)[^/]+$" ) %>% gsub(pattern = ".RData",replacement = "",fixed=T,x=.)
  data_name = str_extract("data/simulation/new_gs/expr_new_gs/expr_rda_rs/expr_ds_1.RData",pattern = "(?<=/)[^/]+$" ) %>% gsub(pattern = ".RData",replacement = "",fixed=T,x=.)
  p_val = get_p_value(data_sce = sce , gwas_zscore = trait_zscore, 
                      gene_mapping_table = mmu_hsa_mapping, group = "cell_ontology_class", 
                      cell_seed_vec = cell_seed_vec)
  return(p_val)
}

#plan(multisession, workers = num_cores)
#result = future_sapply(1:100,sim_p_value_by_idx)

library(parallel)
result = mclapply(1:10000, sim_p_value_by_idx, mc.cores = num_cores)

result_df = data.frame(index = 1:10000, ours_p = unlist(result))
write.table(result_df, file = output_file, sep = "\t",quote = F,row.names = F)
