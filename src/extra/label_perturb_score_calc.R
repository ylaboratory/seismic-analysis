#This script is similar to src/causal_sim.R
#Only save the results of the target cell type
#For each time it will first generate the background 
#1: parameter df file path 
#2: summarized output file header
#3: which column (name) indicates the final output file header?
#4: regular expression contains the pattern of target cell types
#5: number of cores to use

args <- commandArgs(trailingOnly = TRUE)
options(warn = -1)

parameter_df_file <- args[1]
sce_data_path <- args[2]
final_res_path <- args[3]
output_header_col <- args[4]
tar_ct_col <- args[5]
num_cores <- args[6]

#load packages
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("here"))
suppressMessages(library("seismicGWAS"))
suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("parallel"))
suppressMessages(library("limma"))
suppressMessages(library("genefilter"))

#load scripts with functions
source(here("src","tools","sparse_mat_util.R"))
source(here("src","tools","magma_fuma_file_prep.R"))

#load data
load(sce_data_path)
facs_obj_sce <- facs_obj_sce[,!is.na(facs_obj_sce$cell_ontology_class)]
facs_obj_sce$cellid <- paste0("cell_", 1:ncol(facs_obj_sce))
colnames(facs_obj_sce) <- facs_obj_sce$cellid

#calculate logCPM and CPM
assay(facs_obj_sce, "cpm") <- sweep_sparse(assay(facs_obj_sce, "counts")*1e6, margin = 2, stats = colSums(assay(facs_obj_sce, "counts"))) 
assay(facs_obj_sce, "logcpm") <-  assay(facs_obj_sce, "cpm") %>% transform_sparse(function(x) log2(x+1))

#load data and set output path
parameter_df <- read.table(parameter_df_file , header = T, sep = "\t")
process_file <- paste0(final_res_path, "/process.log")

get_seismic_sscore <- function(data_sce, target_gene_list, assay_name = "logcounts", group = "cell_ontology_class", target_cell_type = "target_cell_type"){
  spec_score <- calc_specificity(data_sce, assay_name = assay_name, ct_label_col = group, min_ct_size = 5) 
  gene_idx <- match(target_gene_list, rownames(spec_score)) %>% .[!is.na(.)]
  spec_score <- spec_score %>%
    .[gene_idx, ] %>%
    as.matrix() %>%
    as_tibble(rownames = "gene")
  
  return(spec_score[[target_cell_type]] %>% set_names(spec_score[["gene"]]))
}

#the Bryois et al specificity
get_spc <- function(data_sce, target_gene_list, assay_name = "logcpm", group = "cell_ontology_class", target_cell_type = "target_cell_type"){
  ct_mean <- calc_ct_mean(data_sce, assay_name = "counts", ct_label_col = "cell_ontology_class")
  ct_mean <- sweep(ct_mean*1e6, MARGIN=1, STATS=rowSums(ct_mean), FUN="/")
  spec_score <- sweep(ct_mean, MARGIN = 2, STATS = colSums(ct_mean), FUN = "/")
  spec_score <- t(spec_score[,match(target_gene_list, colnames(spec_score))]) %>%
    as.matrix() %>%
    as_tibble(rownames = "gene")
  
  return(spec_score[[target_cell_type]] %>% set_names(spec_score[["gene"]]))
}

#descore vector
get_de_score <- function(data_sce, target_gene_list, assay_name = "logcpm", group = "cell_ontology_class", target_cell_type = "target_cell_type"){
  #limma fitting
  design <- model.matrix(~factor(data_sce[[group]] ==  target_cell_type))
  colnames(design) <- c("Intercept", target_cell_type)
  fit <- lmFit(assay(data_sce, assay_name), design) %>% eBayes()
  #get results
  de_results <- topTable(fit, coef=target_cell_type, number=Inf) %>%
    as_tibble(rownames = "gene") %>%
    mutate(descore = -log10(adj.P.Val))
  
  de_results <- de_results[match(target_gene_list, de_results[["gene"]]),]
  
  return(de_results[["descore"]] %>% set_names(de_results[["gene"]]))
  
}

#get specificity index
get_si <- function(data_sce, target_gene_list, assay_name = "cpm", group = "cell_ontology_class", target_cell_type = "target_cell_type"){
  ct_mean <- calc_ct_mean(data_sce, assay_name = assay_name, ct_label_col = group) %>%
    .[,match(target_gene_list, colnames(.))] %>%
    t() %>%
    as.matrix()
  
  ct_ratio <- sweep(ct_mean[, colnames(ct_mean) != target_cell_type] + 1e-2, MARGIN = 1, STATS = ct_mean[, target_cell_type] + 1e-2, FUN = "/")
  ct_ratio_rank <- t(matrixStats::colRanks(-ct_ratio, useNames = T, ties.method = "average"))
  
  return(rowMeans(ct_ratio_rank))
  
}



#get all values
get_all_score <- function(data_sce, seed_sample_file, output_header, analysis_group = "cell_ontology_class",  target_cell_type = "target_cell_type"){
  sce <- data_sce
  seed_tbl <- read.table(seed_sample_file, header=T, sep="\t")
  target_gene_list <- rownames(sce)
  
  #calculate score without any modification
  tmp_cell_anno <- colData(sce) %>%
    as.data.frame() %>%
    group_by(sce[[analysis_group]]) %>%
    add_count() %>%
    ungroup() %>%
    filter(n >= 5) 
  
  tmp_sce <- sce[,sce[["cellid"]] %in% tmp_cell_anno[["cellid"]]]
  
  seismic_score <- get_seismic_sscore(tmp_sce, target_gene_list = target_gene_list, target_cell_type = target_cell_type)
  spc_score <- get_spc(tmp_sce, target_gene_list = names(seismic_score), target_cell_type = target_cell_type)
  de_score <- get_de_score(tmp_sce, target_gene_list = names(seismic_score), target_cell_type = target_cell_type)
  si <- get_si(tmp_sce, target_gene_list = names(seismic_score), target_cell_type = target_cell_type)

  res_df <- list(seismic_score, spc_score, de_score, si) %>%
    set_names("seismic_score", "spc_score", "de_score", "si")%>%
    map(~as_tibble(.x, rownames = c("gene"))) %>% 
    map2(names(.), ~set_colnames(.x, c("gene", .y))) %>%
    purrr::reduce(~full_join(.x, .y, by = "gene"))
  

  write.table(res_df, file = paste0(output_header,".all.original_score_df.txt"), row.names = F, col.names = T, sep="\t", quote=F)
  
  
  #extract correct cell anno
  sampled_cell_anno <- colnames(seed_tbl) %>%
    set_names(colnames(seed_tbl)) %>%
    map(~set_names(sce[[analysis_group]], sce$cellid)) %>%
    map2(names(.), ~{tmp = .x; tmp[seed_tbl[[.y]]] = target_cell_type; tmp}) %>%
    map(~as_tibble(.x, rownames = "cellid")) %>%
    map(~group_by(.x, value) %>% add_count() %>% ungroup() %>% filter(n>=5))
  
  all_res <- sampled_cell_anno[1:5] %>%
    map2(names(.)[1:5], ~{
      print(.y)
      tmp_sce <- sce[, sce[["cellid"]] %in% .x[["cellid"]] ];
      tmp_sce[[analysis_group]] <- .x[["value"]][match(tmp_sce[["cellid"]], .x[["cellid"]])];

      seismic_score <- get_seismic_sscore(tmp_sce, target_gene_list = target_gene_list, target_cell_type = target_cell_type)
      spc_score <- get_spc(tmp_sce, target_gene_list = names(seismic_score), target_cell_type = target_cell_type)
      de_score <- get_de_score(tmp_sce, target_gene_list = names(seismic_score), target_cell_type = target_cell_type)
      si <- get_si(tmp_sce, target_gene_list = names(seismic_score), target_cell_type = target_cell_type)

      res_df <- list(seismic_score, spc_score, de_score, si) %>%
        set_names("seismic_score", "spc_score", "de_score", "si")%>%
        map(~as_tibble(.x, rownames = c("gene"))) %>% 
        map2(names(.), ~set_colnames(.x, c("gene", .y))) %>%
        purrr::reduce(~full_join(.x, .y, by = "gene"))
      print("159")
      write.table(res_df, file = paste0(output_header,".all.",.y, ".score_df.txt"), row.names = F, col.names = T, sep="\t", quote=F)
      return(1)
    })
  
  return(1)
}


#function by index
get_value_by_idx <- function(i){
  print(paste0("running ",i))
  tryCatch({
    data_sce <- facs_obj_sce
    output_header <- parameter_df[[output_header_col]][i]
    seed_sample_file <- parameter_df$seed_sample_file[i]
    cell_anno_file <- parameter_df$cell_anno_file[i]
    target_cell_type <- parameter_df[[tar_ct_col]][i]
    result <- get_all_score(data_sce = data_sce, 
                            seed_sample_file = seed_sample_file,
                            output_header = output_header,
                            target_cell_type = target_cell_type)
    return(result)
  }, error = function(e) {
    message(sprintf("Error in processing index %d: %s", i, e$message))
    return(NA)  # Return NA or another appropriate value that signifies an error
  })
}

result <- mclapply(1:nrow(parameter_df), get_value_by_idx, mc.cores = num_cores)
#
print("Finish!")