#load arguments
#1: parameter df file path
#2: output directory path
#3: output directory is organized by which column of the parameter df?
#4: if the influential gene analysis should be done
#5: number of cores used

args <- commandArgs(trailingOnly = TRUE)
options(warn = -1)

parameter_df_file <- args[1]
final_res_path <- args[2]
#organized_by <- args[3]
output_header_col <- args[3]
tar_ct <- args[4]
extract_target_ct <- args[5]
do_inf_analysis <- args[6] #if it was no/false/False/NO, then no inf analysis will be performed. 
#Else if it's a regular expression then do influential analysis for the matched cell type
num_cores <- args[7]


#modify paramters
if(tolower(do_inf_analysis) %in% c("false","f","no","n")){
  do_inf_analysis <- FALSE
}else{
  do_inf_analysis <- TRUE
}

if(tolower(extract_target_ct) %in% c("true","t")){
  extract_target_ct <- TRUE
}else{
  extract_target_ct <- FALSE
}

#load packages
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("here"))
suppressMessages(library("seismicGWAS"))
suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("parallel"))

source(here("src","tools","sparse_mat_util.R"))
source(here("src","tools","magma_fuma_file_prep.R"))

#load data and set output path
parameter_df <- read.table(parameter_df_file , header = T)
process_file <- paste0(final_res_path, "/process.log")
final_output_file <- paste0(final_res_path, "/all_res.txt")
mmu_hsa_mapping_file <- here("data","ref","mapping","mmu_hsa_mapping.rda")

load(mmu_hsa_mapping_file)  
mmu_hsa_mapping <- mmu_hsa_mapping %>% 
  distinct(mmu_symbol, hsa_entrez) %>%
  drop_na() %>%
  group_by(mmu_symbol) %>% 
  filter(n()==1) %>%
  group_by(hsa_entrez) %>%
  filter(n()==1) %>%
  ungroup() %>%
  mutate(hsa_entrez = as.character(hsa_entrez))


perturb_sce <- function(data_sce, replace_mat, new_cell_anno, group = "cell_ontology_class"){
  #change cell type annotation
  colData(data_sce)[[group]][match(names(new_cell_anno), colnames(data_sce))] = new_cell_anno
  
  #get matched indices
  row_indices <- match(rownames(replace_mat), rownames(data_sce))
  col_indices <- match(colnames(replace_mat), colnames(data_sce))
  
  #modify count matrix
  counts_mat <- as.matrix(assay(data_sce, "counts"))
  size_factor = colSums(counts_mat) / colSums(2^assay(data_sce, "logcounts") -1)  #size factor
  counts_mat[row_indices, col_indices] <- replace_mat
  counts_mat <- as(counts_mat, "dgCMatrix")
  
  #other normalization
  logcounts_mat <- sweep_sparse(counts_mat, margin = 2, stats = size_factor) %>% transform_sparse(function(x) log2(x+1))
  logcpm_mat <- sweep_sparse(counts_mat*1e6, margin = 2, stats = colSums(counts_mat)) %>% transform_sparse(function(x) log2(x+1))
  
  assay(data_sce, "counts") <- counts_mat
  assay(data_sce, "logcounts") <- logcounts_mat
  assay(data_sce, "logcpm") <- logcpm_mat
  
  return(data_sce)
}


#function for seismic analysis
seismic_fdr <- function(data_sce, assay_name = "logcounts", gwas_zscore_df, output_header, group = "cell_ontology_class", target_cell_type_list = "target_cell_type", influential_analysis = T){
  #seismic workflow
  sscore <- calc_specificity(data_sce, assay_name = assay_name, ct_label_col = group)
  sscore_hsa <- translate_gene_ids(sscore, from = "mmu_symbol")
  all_association <- get_ct_trait_associations(sscore = sscore_hsa, magma = gwas_zscore_df)
  
  #target_cell_type_idx 
  target_ct_idx = match(target_cell_type_list, all_association$cell_type)
  
  #influential analysis
  if(influential_analysis){
    for (cur_ct_idx in target_ct_idx[which(!is.na(target_ct_idx))]){
      example_inf_genes <- find_inf_genes(all_association$cell_type[cur_ct_idx], sscore_hsa, magma = gwas_zscore_df)
      write.table(example_inf_genes, file = paste0(output_header,".", all_association$cell_type[cur_ct_idx] ,".inf_genes.txt"), col.names = T, row.names = F, quote = F, sep="\t")
    }
  }
  
  #save results
  write.table(all_association, file = paste0(output_header, ".results.txt"), col.names = T, row.names = F, quote = F, sep="\t")
  
  if (!extract_target_ct){
    return(rep(NA, length(target_cell_type_list)))
  }
  
  #return
  return(all_association$FDR[target_ct_idx])
}

#function for S-MAGMA analysis
magma_fdr <- function(data_sce, assay_name = "counts", magma_raw_path, output_header, gene_mapping_tbl, group = "cell_ontology_class", target_cell_type_list = "target_cell_type"){
  mean_mat <- calc_ct_mean(data_sce, assay_name = assay_name, ct_label_col = group)
  
  mean_mat = mean_mat[, colnames(mean_mat) %in% gene_mapping_tbl$mmu_symbol] %>% 
    set_colnames(gene_mapping_tbl$hsa_entrez[match( colnames(.), gene_mapping_tbl$mmu_symbol)])
  
  mean_mat = mean_mat[, which(colSums(mean_mat)>0)]
  
  mean_mat <- sweep(mean_mat*1e6, MARGIN=1, STATS=rowSums(mean_mat), FUN="/")
  
  print_magma_fuma_tbl(mean_mat, "MAGMA", 
                       main_table_path = paste0(output_header, ".magma.txt"),
                       aux_table_path = paste0(output_header, ".magma.aux.txt"), verbose = F)
  
  #run magma
  ms <- system(paste0(here("bin", "magma", "magma"), " --gene-results ", magma_raw_path, " --set-annot ", output_header, ".magma.txt --out ", output_header,".magma"), intern = TRUE)
  
  #extract results or not
  if (!extract_target_ct){
    return(rep(NA, length(target_cell_type_list)))
  }
  
  #read in results
  magma_res_tbl <- read.table(paste0(output_header, ".magma.gsa.out"), header = TRUE)
  magma_annot_tbl <- read.table(paste0(output_header, ".magma.aux.txt"), header = TRUE, sep="\t")
  
  #get the encoded name of the cell type
  encoded_ct_list <- magma_annot_tbl$encoded_name[match(target_cell_type_list, magma_annot_tbl$cell_type )] 
  
  #get results
  magma_res_tbl <- magma_res_tbl %>%
    mutate(fdr = p.adjust(P, method="fdr"))
  res_fdr <- magma_res_tbl$fdr[match(encoded_ct_list, magma_res_tbl$VARIABLE)]
  
  return(res_fdr)
}

#function for FUMA analysis
fuma_fdr <- function(data_sce, assay_name="logcpm", magma_raw_path, output_header, group = "cell_ontology_class", target_cell_type_list = "target_cell_type"){
  #prepare fuma file
  mean_mat <- calc_ct_mean(data_sce, assay_name =  assay_name, ct_label_col = group)
  mean_mat_hsa <- translate_gene_ids(t(mean_mat), from = "mmu_symbol")
  print_magma_fuma_tbl(t(mean_mat_hsa), "FUMA", 
                       main_table_path = paste0(output_header, ".fuma.txt"),
                       aux_table_path = paste0(output_header, ".fuma.aux.txt"), verbose = F)
  
  ms = system(paste0(here("bin", "magma", "magma"), " --gene-results ", magma_raw_path, " --gene-covar ", output_header, ".fuma.txt --model condition-hide=Average direction=greater --out ", output_header,".fuma"), intern = TRUE)
  
  #extract results or not
  if (!extract_target_ct | !file.exists(paste0(output_header, ".fuma.gsa.out"))){ #if the result does not exist due to 
    return(rep(NA, length(target_cell_type_list)))
  }
  
  #read in results
  fuma_res_tbl <- read.table(paste0(output_header, ".fuma.gsa.out"), header = TRUE)
  fuma_annot_tbl <- read.table(paste0(output_header, ".fuma.aux.txt"), header = TRUE, sep="\t")
  
  #get the encoded name of the cell type
  encoded_ct_list <- fuma_annot_tbl$encoded_name[match(target_cell_type_list, fuma_annot_tbl$cell_type )] 
  
  #get results
  fuma_res_tbl <- fuma_res_tbl %>%
    mutate(fdr = p.adjust(P, method="fdr"))
  res_fdr <- fuma_res_tbl$fdr[match(encoded_ct_list, fuma_res_tbl$VARIABLE)]
  
  return(res_fdr)
}

#return a vector of fdr
get_fdr <- function(data_sce_file, replace_mat_file, cell_anno_file, zscore_file, magma_raw_file,  output_header, group_to_replace = "new_cell_ontology_class", gene_mapping_table = mmu_hsa_mapping, target_cell_type_pattern = "target_cell_type", influential_analysis = T) {
  load(data_sce_file)
  replace_mat <- read.table(replace_mat_file, header = T) %>% as.matrix
  cell_anno <- read.table(cell_anno_file, header = T, sep="\t")
  gwas_zscore_df <- read.table(zscore_file, header = T)
  
  #perturbation
  colnames(sce) <- sce$cellid
  cellid_idx <- match(colnames(replace_mat), cell_anno$cellid)
  new_cell_anno <- cell_anno[[group_to_replace]][cellid_idx] %>% set_names(cell_anno$cellid[cellid_idx])
  sce <- perturb_sce(sce, replace_mat = replace_mat, new_cell_anno = new_cell_anno)
  
  #extract target cell type names by regular expression
  target_cell_type_list <- unique(cell_anno[[group_to_replace]])
  target_cell_type_list <-target_cell_type_list[grepl(pattern =  target_cell_type_pattern, x=target_cell_type_list)]
  
  seismic_fdr_value <- seismic_fdr(sce, gwas_zscore_df = gwas_zscore_df, output_header = output_header, influential_analysis = influential_analysis, target_cell_type_list = target_cell_type_list)
  magma_fdr_value <- magma_fdr(sce,  magma_raw_path = magma_raw_file, output_header = output_header, gene_mapping_tbl = gene_mapping_table, target_cell_type_list = target_cell_type_list)
  fuma_fdr_value <- fuma_fdr(sce, magma_raw_path = magma_raw_file, output_header = output_header, target_cell_type_list = target_cell_type_list)
  
  #when FUMA return nothing because of correlation
  if (is.null(fuma_fdr_value )) {
    fuma_fdr_value = rep(NA, length(target_cell_type_list))
  }

  return(list(target_cell_type_list, c(seismic_fdr_value, magma_fdr_value, fuma_fdr_value)))
}

get_fdr_by_idx <- function(i){
  if (i %% 10 == 0) {
    cat(paste("Running seed ", i, "\n"), file = process_file, append = TRUE)
  }
  tryCatch({
    data_sce_file <- parameter_df$expr_rda_file[i]
    replace_mat_file <- parameter_df$perturbed_mat_file[i]
    cell_anno_file <- parameter_df$cell_anno_file[i]
    zscore_file <- parameter_df$gs_zscore_file[i]
    magma_raw_file <- parameter_df$magma_raw_file[i]
    output_header <- parameter_df[[output_header_col]][i]
    result <- get_fdr(data_sce_file = data_sce_file, 
                      replace_mat_file = replace_mat_file, 
                      cell_anno_file = cell_anno_file,
                      zscore_file = zscore_file, 
                      magma_raw_file = magma_raw_file, 
                      output_header = output_header,
                      target_cell_type_pattern = tar_ct,
                      influential_analysis = do_inf_analysis)
    return(result)
  }, error = function(e) {
    message(sprintf("Error in processing index %d: %s", i, e$message))
    return(NA)  # Return NA or another appropriate value that signifies an error
  })
}


result <- mclapply(1:nrow(parameter_df), get_fdr_by_idx, mc.cores = num_cores)
#result <- mclapply(1:50, get_fdr_by_idx, mc.cores = num_cores)

result_df <- data.frame(index = 1:nrow(parameter_df), purrr::list_transpose(map(result, ~.x[[2]]))) %>% 
   set_names(c("index", paste0(result[[1]][[1]], ".seismic_fdr"),  paste0(result[[1]][[1]], ".magma_fdr"), paste0(result[[1]][[1]], ".fuma_fdr")))
# result_df <- data.frame(index = 1:50, purrr::list_transpose(map(result, ~.x[[2]]))) %>% 
#   set_names(c("index", paste0(result[[1]][[1]], ".seismic_fdr"),  paste0(result[[1]][[1]], ".magma_fdr"), paste0(result[[1]][[1]], ".fuma_fdr")))


write.table(result_df, file = final_output_file, sep = "\t",quote = F, row.names = F)
