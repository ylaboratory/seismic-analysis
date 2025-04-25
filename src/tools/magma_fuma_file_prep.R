#function to calculate mean: (truncated from calc_specificity)
calc_ct_mean <- function(sce, assay_name = "logcounts",
                             ct_label_col = "idents",
                             min_ct_size = 20) {
  ct <- N <- nz.count <- ave_exp_ct <- NULL # due to non-standard evaluation notes in R CMD check
  
  # data formatting checks
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("Only SingleCellExperiment class input is accepted. Please supply an
         SingleCellExperiment object and try again.")
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(sce)) {
    stop("Assay '", assay_name, "' does not exist in the
         SingleCellExperiment object. Please choose a valid
         assay in the SingleCellExperiment.")
  }
  
  if (!ct_label_col %in% names(SummarizedExperiment::colData(sce))) {
    stop("ct_label_col '", ct_label_col, "' does not exist in the
         SingleCellExperiment object. Please choose a valid
         column of cell type labels in the SingleCellExperiment.")
  }
  
  # extract the log normalized counts (dgCMatrix)
  data_mat <- SummarizedExperiment::assay(sce, assay_name)
  
  # extract the cell metadata
  cell_meta <- SummarizedExperiment::colData(sce)
  
  # make sure the cell name exist
  if (any(is.null(rownames(cell_meta))) || any(is.null(colnames(data_mat)))){
    rownames(cell_meta) <- colnames(data_mat) <- paste0("cell.",1:ncol(data_mat))
  }
  
  # extract cell type grouping
  ct_groups <- data.table::data.table(
    cell = rownames(cell_meta),
    ct = cell_meta[[ct_label_col]], key = "ct"
  )
  
  # check that there are at least a few different cell types
  ct_groups_n <- ct_groups[, .N, by = ct]
  
  # filter out cell types that do not have minimum # of cells
  ct_groups_n <- ct_groups_n[N >= min_ct_size]
  ct_groups <- ct_groups[ct %in% ct_groups_n$ct]
  data_mat <- data_mat[, ct_groups$cell]
  
  # filter out genes that do not have minimal cell coverage
  stats.dt <- data.table::as.data.table(Matrix::rowSums(data_mat != 0), keep.rownames = T) %>%
    magrittr::set_colnames(c("gene", "nz.count"))
  
  data_mat <- data_mat[stats.dt$gene, ]
  
  # calculate mean gene expression per cell type
  factor_mat <- Matrix::fac2sparse(factor(ct_groups$ct, levels = unique(ct_groups$ct)))
  sum_mat <- Matrix::t(data_mat %*% Matrix::t(factor_mat))
  #return(list(factor_mat, data_mat, ct_groups))
  mean_mat <- sum_mat %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat)) %>%
    magrittr::set_rownames(unique(ct_groups$ct))
  
  
  return(as(mean_mat, "unpackedMatrix"))
}


#function to calculate mean: (truncated from calc_specificity)
calc_ct_ratio <- function(sce, assay_name = "logcounts",
                         ct_label_col = "idents",
                         min_ct_size = 20) {
  ct <- N <- nz.count <- ave_exp_ct <- NULL # due to non-standard evaluation notes in R CMD check
  
  # data formatting checks
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("Only SingleCellExperiment class input is accepted. Please supply an
         SingleCellExperiment object and try again.")
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(sce)) {
    stop("Assay '", assay_name, "' does not exist in the
         SingleCellExperiment object. Please choose a valid
         assay in the SingleCellExperiment.")
  }
  
  if (!ct_label_col %in% names(SummarizedExperiment::colData(sce))) {
    stop("ct_label_col '", ct_label_col, "' does not exist in the
         SingleCellExperiment object. Please choose a valid
         column of cell type labels in the SingleCellExperiment.")
  }
  
  # extract the log normalized counts (dgCMatrix)
  data_mat <- SummarizedExperiment::assay(sce, assay_name)
  
  # extract the cell metadata
  cell_meta <- SummarizedExperiment::colData(sce)
  
  # make sure the cell name exist
  if (any(is.null(rownames(cell_meta))) || any(is.null(colnames(data_mat)))){
    rownames(cell_meta) <- colnames(data_mat) <- paste0("cell.",1:ncol(data_mat))
  }
  
  # extract cell type grouping
  ct_groups <- data.table::data.table(
    cell = rownames(cell_meta),
    ct = cell_meta[[ct_label_col]], key = "ct"
  )
  
  # check that there are at least a few different cell types
  ct_groups_n <- ct_groups[, .N, by = ct]
  
  # filter out cell types that do not have minimum # of cells
  ct_groups_n <- ct_groups_n[N >= min_ct_size]
  ct_groups <- ct_groups[ct %in% ct_groups_n$ct]
  data_mat <- data_mat[, ct_groups$cell]
  
  # filter out genes that do not have minimal cell coverage
  stats.dt <- data.table::as.data.table(Matrix::rowSums(data_mat != 0), keep.rownames = T) %>%
    magrittr::set_colnames(c("gene", "nz.count"))
  
  data_mat <- data_mat[stats.dt$gene, ]
  
  # calculate mean gene expression per cell type
  factor_mat <- Matrix::fac2sparse(factor(ct_groups$ct, levels = unique(ct_groups$ct)))
  
  # calculate ratio of cells a gene has non-zero expression in per cell type
  ratio_mat <- Matrix::t((data_mat > 0) %*% Matrix::t(factor_mat)) %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat)) %>%
    magrittr::set_rownames(unique(ct_groups$ct))
  
  return(as(ratio_mat, "unpackedMatrix"))
}

#function to print out magma and fuma tables for later analysis
print_magma_fuma_tbl <- function(mean_mat, table_type, main_table_path,  aux_table_path = NULL, verbose = T) {
  #main table directory
  if (!dir.exists(base::dirname(main_table_path))) {
    message(paste0("The directory of the main table path:", base::dirname(main_table_path) ," does not exist... cerated one.!"))
    dir.create(base::dirname(main_table_path), recursive = TRUE)
  }
  
  #if output MAGMA table
  if(table_type == "MAGMA"){
    main_tbl <- sweep(mean_mat, MARGIN = 2, STATS = Matrix::colSums(mean_mat), FUN="/") %>% 
      as.matrix() %>%
      Matrix::t() %>% 
      dplyr::as_tibble(rownames = "hsa_entrez") %>%
      tidyr::pivot_longer(!hsa_entrez, names_to = "cell_type", values_to = "specificity") %>%
      dplyr::group_by(cell_type) %>%
      dplyr::slice_max(specificity, prop=0.1) %>% 
      dplyr::summarize( genes = paste(hsa_entrez, collapse = " ")) %>%
      dplyr::arrange(match(cell_type, rownames(mean_mat))) %>%
      dplyr::mutate(cell_type = paste0("cluster.",1:n()))
    write.table(main_tbl,file=main_table_path, col.names = F, row.names = F, sep=" ",quote=F)
  }else{
    #if output FUMA table
    main_tbl <- rbind(mean_mat, Matrix::colMeans(mean_mat)) %>%
      magrittr::set_rownames(c(rownames(mean_mat), "Average")) %>%
      as.matrix() %>%
      Matrix::t() %>% 
      dplyr::as_tibble(rownames = "GENE") %>%
      magrittr::set_colnames(c("GENE", paste0("cluster.",1:(ncol(.)-2) ), "Average" ))
    write.table(main_tbl,file=main_table_path, col.names =  T, row.names = F, sep=" ",quote=F)
  }
  if(verbose){
    message(paste0("The main table has been printed out as ",main_table_path))
  }
  
  
  #if print auxiliary file
  if(!is.null(aux_table_path)){
    aux_tbl <- dplyr::tibble(cell_type = rownames(mean_mat)) %>%
      dplyr::mutate(encoded_name = paste0("cluster.",1:n())) 
    #print the table
    if (!dir.exists(base::dirname(aux_table_path))) {
      message(paste0("The directory of the auxilary table path:", base::dirname(aux_table_path) ," does not exist... cerated one.!"))
      dir.create(base::dirname(aux_table_path), recursive = TRUE)
    }
    write.table(aux_tbl,file=aux_table_path, col.names =  T, row.names = F, sep="\t",quote=F)
    if(verbose){
      message(paste0("The auxiliary table has been printed out as ",aux_table_path))
    }
  }
}
