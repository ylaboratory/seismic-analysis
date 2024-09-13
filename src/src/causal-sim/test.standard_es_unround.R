#causal simulation:
#load sce, perturbed matrix and print out results
#1: sce directory
#2: sce directory
#3: MAGMA Z-score file_directory 
#4: MAGMA raw file directory
#5: Final output file directory
#6: Number of cores


#perturbed expression matrix are encoded as es"$EFFECT_SIZE"_gs"$GENE_SET_NUM"_ds"$DATASET_NUM"
#the script will try to match gene set number and dataset number in the directory
#load arguments

#load packages
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("here"))
suppressMessages(library("seismicGWAS"))
suppressMessages(library("SingleCellExperiment"))

#load source code
source(here("src","tools","sparse_mat_util.R"))
source(here("src","causal-sim","causal_sim_func.R"))


#load data
parameter_df <- read.table(here("data", "expr", "causal_sim", "standard", "perturbed_expr_unround", "parameter_df.txt"), header = T)
output_dir <- here("results","causal_sim","standard_unround")
mmu_hsa_mapping <- here("data","ref","mapping","mmu_hsa_mapping.rda")
process_file <- paste0(output_dir, "/process.log")
num_cores <- 15

#load mapping
load(mmu_hsa_mapping)
mmu_hsa_mapping <- mmu_hsa_mapping %>% 
  distinct(mmu_symbol, hsa_entrez) %>%
  drop_na() %>%
  group_by(mmu_symbol) %>% 
  filter(n()==1) %>%
  group_by(hsa_entrez) %>%
  filter(n()==1) %>%
  ungroup() %>%
  mutate(hsa_entrez = as.character(hsa_entrez))

#return a vector of fdr
get_fdr <- function(data_sce_file, replace_mat_file, zscore_file, magma_raw_file,  output_header, gene_mapping_table = mmu_hsa_mapping) {
  load(data_sce_file)
  replace_mat <- read.table(replace_mat_file, header = T) %>% as.matrix
  gwas_zscore_df <- read.table(zscore_file, header = T)
  
  colnames(sce) <- sce$cellid
  sce <- perturb_sce(sce, replace_mat = replace_mat)
  
  seismic_fdr_value <- seismic_fdr(sce, gwas_zscore_df = gwas_zscore_df, output_header = output_header)
  magma_fdr_value <- magma_fdr(sce,  magma_raw_path = magma_raw_file, output_header = output_header, gene_mapping_tbl = gene_mapping_table)
  fuma_fdr_value <- fuma_fdr(sce, magma_raw_path = magma_raw_file, output_header = output_header)
  
  #when FUMA return nothing because of correlation
  if (is.null(fuma_fdr_value )) {
    fuma_fdr_value = NA
  }
  
  #return
  return(c(seismic_fdr_value, magma_fdr_value, fuma_fdr_value))
}

get_fdr_by_idx <- function(i){
  if (i %% 10 == 0) {
    cat(paste("Running seed ", i, "\n"), file = process_file, append = TRUE)
  }
  tryCatch({
    data_sce_file <- parameter_df$expr_rda_file[i]
    replace_mat_file <- parameter_df$perturbed_mat_file[i]
    zscore_file <- parameter_df$gs_zscore_file[i]
    magma_raw_file <- parameter_df$magma_raw_file[i]
    output_header = paste0(output_dir, "/",parameter_df$es_name[i],"/",parameter_df$gs_name[i], ".", parameter_df$expr_name[i])
    result <- get_fdr(data_sce_file = data_sce_file, 
                      replace_mat_file = replace_mat_file, 
                      zscore_file = zscore_file, 
                      magma_raw_file = magma_raw_file, 
                      output_header = output_header)
    return(result)
  }, error = function(e) {
    message(sprintf("Error in processing index %d: %s", i, e$message))
    return(NA)  # Return NA or another appropriate value that signifies an error
  })
}

library(parallel)
result <- mclapply(1:2000, get_fdr_by_idx, mc.cores = num_cores)

result_df <- data.frame(index = 1:2000, expr_name = parameter_df$expr_name[1:2000], gs_name = parameter_df$gs_name[1:2000], es_name = parameter_df$es_name[1:2000],
                        seismic_fdr = unlist(map(result, ~.x[1])),magma_fdr = unlist(map(result, ~.x[2])),fuma_fdr = unlist(map(result, ~.x[3])))

write.table(result_df, file = paste0(output_dir, "/all_res.txt"), sep = "\t",quote = F, row.names = F)