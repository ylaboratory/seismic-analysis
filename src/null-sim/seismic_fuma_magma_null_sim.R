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

#load source code
source(here("src","tools","magma_fuma_file_prep.R"))
source(here("src","tools","sparse_mat_util.R"))

#load data
load(args[1])
cell_seed = read.table(args[2], header=T,sep="\t")
z_score_file = args[3]
magma_raw_file = args[4]
process_file = args[5]
num_cores= as.numeric(args[6])
output_file = args[7]
temp_dir = args[8]

load("data/expr/null_sim/expr_rda_rs/expr_ds_1.rda")
cell_seed = read.table("data/expr/null_sim/seed_table/sample_cell_idx.1.txt", header=T,sep="\t")
z_score_file = "data/gwas/null_sim/zscore/gs_1.genes.out"
magma_raw_file = "data/gwas/null_sim/magma_raw/gs_1.genes.raw"
process_file = "data/test.txt"
num_cores = 10
output_file = "data/test_res.txt"
temp_dir = "data/temp"
#load gwas zscore and gene id mapping
trait_zscore = load_zscore(z_score_file)
data("mmu_hsa_mapping")  

#return a data_sce
seismic_p_value = function(data_sce, gwas_zscore, gene_mapping_table, group) {
  data_sce = cal_stat(data_obj = data_sce, meta_data = as.data.frame(colData(data_sce)), group = group, assay_name = "logcounts") #this is log normalized counts 
  data_sce = cal_sscore(data_obj = data_sce) 
  data_sce = trans_mmu_to_hsa_stat(data_sce, gene_mapping_table = gene_mapping_table, from = "mmu_symbol", to = "hsa_entrez")

  data_sce = add_glob_stats(data_sce, stats = c("det_cell_num", "ave_exp_ct", "max_exp_ct")) 
  data_sce = cal_ct_asso(data_sce, gwas_zscore, gene_filter_setting = "det_cell_num>=10 & ave_exp_ct > 0.1 & max_exp_ct > 0.1")
  p_value = get_ct_asso(data_sce, trait_name = "all", asso_model = "linear")[[1]] %>% filter(cell_type == "fake_cell_type") %>% pull(Pvalue)
  
  return(p_value)
}


magma_p_value = function(data_sce, gene_mapping_table, temp_file_header, magma_raw_path, group) {
  data_sce = cal_stat(data_obj = data_sce, meta_data = as.data.frame(colData(data_sce)), group = group, assay_name = "cpm",mean_only = T)
  data_sce = trans_mmu_to_hsa_stat(data_sce, gene_mapping_table = gene_mapping_table, from = "mmu_symbol", to = "hsa_entrez")
  
  # magma_mean_tbl = t(metadata(data_sce)[["seismicGWAS.data"]][["group_info"]][["mean_mat"]]) %>% 
  #   as.matrix() %>%
  #   as_tibble(rownames = "hsa_entrez") %>%
  #   pivot_longer(!hsa_entrez, names_to = "cluster_name", values_to = "specificity") %>%
  #   group_by(cluster_name) %>%
  #   slice_max(specificity, prop = 0.1) %>% 
  #   summarize(genes = paste(hsa_entrez, collapse = " "))
  print_magma_fuma_tbl(data_sce, table_type = "MAGMA", main_table_path =  paste0(temp_file_header, ".magma.txt"), aux_table_path = paste0(temp_file_header, ".magma.aux.txt"))
  
  # write.table(magma_mean_tbl %>% filter(cluster_name == "fake_cell_type"), file = paste0(temp_file_header, ".magma.txt"), col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
  
  ms = system(paste0(here("bin", "magma", "magma"), " --gene-results ", magma_raw_path, " --set-annot ", temp_file_header, ".magma.txt --out ", temp_file_header,".magma"), intern = TRUE)
  
  magma_res_tbl = read.table(paste0(temp_file_header, ".magma.gsa.out"), header = TRUE)
  magma_annot_tbl = read.table(paste0(temp_file_header, ".magma.aux.txt"), header = TRUE, sep="\t")
  
  target_cell_type = magma_annot_tbl %>% 
    filter(cell_type=="fake_cell_type") %>% 
    pull(encoded_name)
  
  ms = system(paste0("rm ", temp_file_header, ".*"), intern = TRUE)
  
  res_p_value = magma_res_tbl %>%
    filter(VARIABLE == target_cell_type) %>%
    pull(P)
  
  return(res_p_value)
}

fuma_p_value = function(data_sce, gene_mapping_table, temp_file_header, magma_raw_path, group) {
  data_sce = cal_stat(data_obj = data_sce, meta_data = as.data.frame(colData(data_sce)), group = group, assay_name = "logcpm", mean_only = T)
  data_sce = trans_mmu_to_hsa_stat(data_sce, gene_mapping_table = gene_mapping_table, from = "mmu_symbol", to = "hsa_entrez")
  
  print_magma_fuma_tbl(data_sce, table_type = "FUMA", main_table_path =  paste0(temp_file_header, ".fuma.txt"), aux_table_path = paste0(temp_file_header, ".fuma.aux.txt"))
  
  ms = system(paste0(here("bin", "magma", "magma"), " --gene-results ", magma_raw_path, " --gene-covar ", temp_file_header, ".fuma.txt --model condition-hide=Average direction=greater --out ", temp_file_header,".fuma"), intern = TRUE)
  
  fuma_res_tbl = read.table(paste0(temp_file_header, ".fuma.gsa.out"), header = TRUE)
  fuma_annot_tbl = read.table(paste0(temp_file_header, ".fuma.aux.txt"), header = TRUE, sep="\t")
  
  target_cell_type = fuma_annot_tbl %>% 
    filter(cell_type=="fake_cell_type") %>% 
    pull(encoded_name)
  
  ms = system(paste0("rm ", temp_file_header, ".*"), intern = TRUE)
  
  res_p_value = fuma_res_tbl %>%
    filter(VARIABLE == target_cell_type) %>%
    pull(P)

  
  return(res_p_value)
}

#return a vector of p value
get_p_value = function(data_sce, gwas_zscore, magma_raw_path, gene_mapping_table, group, cell_seed_vec, temp_file_header) {
  data_sce = reset_seismic_analysis(data_sce)
  colData(data_sce)[[group]][cell_seed_vec] = "fake_cell_type" #fake cell type
  seismic_p = seismic_p_value(data_sce, gwas_zscore = gwas_zscore, gene_mapping_table = gene_mapping_table, group = group)
  magma_p = magma_p_value(data_sce, gene_mapping_table, temp_file_header = temp_file_header, magma_raw_path = magma_raw_path, group = group)
  fuma_p = fuma_p_value(data_sce, gene_mapping_table, temp_file_header = temp_file_header, magma_raw_path = magma_raw_path, group = group)
  
  if (is.null(fuma_p)) {
    fuma_p = NA
  }
  
  return(c(seismic_p, magma_p, fuma_p))
}

#function wrapper
sim_p_value_by_idx = function(i) {
  if (i < 1 | i > nrow(cell_seed)) {
    stop("i is not correct")
  }
  
  if (i %% 10 == 0) {
    cat(paste("Running seed ", i, "\n"), file = process_file, append = TRUE)
  }
  
  cell_seed_vec = cell_seed[i, ] %>% unlist
  data_name = str_extract(args[1], pattern = "(?<=/)[^/]+$") %>% 
    gsub(pattern = ".rda", replacement = "", fixed = TRUE, x = .) #get the rda name 
  
  p_val_vec = get_p_value(data_sce = sce, 
                          gwas_zscore = trait_zscore, 
                          gene_mapping_table = mmu_hsa_mapping,
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
