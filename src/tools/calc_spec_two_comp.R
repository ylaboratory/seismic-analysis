
calc_spec_two_comp <- function(sce, assay_name = "logcounts",
                             ct_label_col = "idents", min_uniq_ct = 2,
                             min_ct_size = 20, min_cells_gene_exp = 10,
                             min_avg_exp_ct = 0.1) {
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
  if (any(is.null(rownames(cell_meta))) || any(is.null(colnames(data_mat)))) {
    rownames(cell_meta) <- colnames(data_mat) <- paste0("cell.",1:ncol(data_mat))
  }
  
  # extract cell type grouping
  ct_groups <- data.table::data.table(
    cell = rownames(cell_meta),
    ct = cell_meta[[ct_label_col]], key = "ct"
  )
  
  # check that there are at least a few different cell types
  ct_groups_n <- ct_groups[, .N, by = ct]
  
  if (nrow(ct_groups_n) < min_uniq_ct) {
    stop("There are fewer than ", min_uniq_ct, " in the SingleCellExperiment.
         Decrease the min_uniq_ct threshold or select cell type labels with
         more unique cell types.")
  }
  
  # filter out cell types that do not have minimum # of cells
  ct_groups_n <- ct_groups_n[N >= min_ct_size]
  ct_groups <- ct_groups[ct %in% ct_groups_n$ct]
  data_mat <- data_mat[, ct_groups$cell]
  
  if (nrow(ct_groups_n) < min_uniq_ct) {
    stop("There are fewer than ", min_uniq_ct, " in the SingleCellExperiment after
          filtering out cell types with fewer than ", min_ct_size, " cells.
         Decrease the min_uniq_ct or min_ct_size thresholds or select cell type labels with
         more unique cell types / larger numbers of cells.")
  }
  
  # filter out genes that do not have minimal cell coverage
  stats.dt <- data.table::as.data.table(Matrix::rowSums(data_mat != 0), keep.rownames = T) %>%
    magrittr::set_colnames(c("gene", "nz.count"))
  stats.dt <- stats.dt[nz.count >= min_cells_gene_exp]
  
  if (nrow(stats.dt) < 1000) {
    warning("There are fewer than 1000 genes in the SingleCellExperiment that
            are expressed in at least ", min_cells_gene_exp, " cells.
            Consider relaxing the threshold or double check the input file.")
  }
  
  data_mat <- data_mat[stats.dt$gene, ]
  
  # calculate mean gene expression per cell type
  factor_mat <- Matrix::fac2sparse(factor(ct_groups$ct, levels = unique(ct_groups$ct)))
  sum_mat <- Matrix::t(data_mat %*% Matrix::t(factor_mat))
  mean_mat <- sum_mat %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat)) %>%
    magrittr::set_rownames(unique(ct_groups$ct))
  
  # filter out genes whose cell-type averaged expression does not exceed a baseline level
  stats.dt$ave_exp_ct <- Matrix::colMeans(mean_mat)
  stats.dt <- stats.dt[ave_exp_ct >= min_avg_exp_ct]
  
  if (nrow(stats.dt) < 1000) {
    warning("There are fewer than 1000 genes in the SingleCellExperiment that
            do not have at least mean expression ", min_avg_exp_ct, " across cell types.
            Consider relaxing the threshold or double check the input file.")
  }
  
  data_mat <- data_mat[stats.dt$gene, ]
  sum_mat <- sum_mat[, stats.dt$gene]
  mean_mat <- mean_mat[, stats.dt$gene]
  
  # calculate variance of gene expression per cell type
  var_mat <- (Matrix::t(data_mat^2 %*% Matrix::t(factor_mat)) - 2 * mean_mat * sum_mat +
                sweep_sparse(mean_mat^2, margin = 1, stats = ct_groups_n$N, fun = "*")) %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N - 1, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat)) %>%
    magrittr::set_rownames(unique(ct_groups$ct))
  
  # calculate ratio of cells a gene has non-zero expression in per cell type
  ratio_mat <- Matrix::t((data_mat > 0) %*% Matrix::t(factor_mat)) %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat)) %>%
    magrittr::set_rownames(unique(ct_groups$ct))
  
  # calculate an indicator matrix for all other cell types
  out_group_mat <- matrix(1, nrow = dim(ct_groups_n)[1], ncol = dim(ct_groups_n)[1]) -
    diag(nrow = dim(ct_groups_n)[1], ncol = dim(ct_groups_n)[1])
  
  # mean for out group per cell type
  out_mean <- (out_group_mat %*% sum_mat) %>%
    sweep_sparse(margin = 1, stats = out_group_mat %*% ct_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat)) %>%
    magrittr::set_rownames(unique(ct_groups$ct))
  
  # variance for out groups per cell type
  tot_var <- sweep_sparse(var_mat, margin = 1, stats = ct_groups_n$N - 1, fun = "*")
  tot_mean_sq <- sweep_sparse(mean_mat^2, margin = 1, stats = ct_groups_n$N, fun = "*")
  out_variance <- (out_group_mat %*% tot_var + out_group_mat %*% tot_mean_sq -
                     sweep_sparse(out_mean^2, margin = 1, stats = out_group_mat %*% ct_groups_n$N, fun = "*")) %>%
    sweep_sparse(margin = 1, stats = out_group_mat %*% ct_groups_n$N - 1, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat)) %>%
    magrittr::set_rownames(unique(ct_groups$ct))
  
  # probability of relatively higher expression for each gene in each cell type (vs other genes)
  rel_exp <- (mean_mat - out_mean) /
    sqrt(sweep_sparse(var_mat, margin = 1, stats = ct_groups_n$N - 1, fun = "/") +
           sweep_sparse(out_variance, margin = 1, stats = out_group_mat %*% ct_groups_n$N - 1, fun = "/"))
  
  # calculate seismic specificity score
  spec_score <- stats::pnorm(as.matrix(rel_exp)) * ratio_mat
  spec_score <- sweep_sparse(x = spec_score, margin = 2, stats = Matrix::colSums(spec_score), fun = "/")
  
  return(list(spec_score = Matrix::t(spec_score),
              ratio = Matrix::t(ratio_mat),
              p = t(stats::pnorm(as.matrix(rel_exp)))))}