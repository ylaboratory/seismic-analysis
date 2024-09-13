library(seismicGWAS)
parameter_df <- read.table("data/expr/causal_sim/sc3_gene_num/gene_400/perturbed_expr/parameter_df.txt", header = T) %>% as_tibble()
i <- 1

#load data
load(parameter_df$expr_rda_file[i])
data_sce <- sce
colnames(data_sce) <- data_sce$cellid

cell_anno <- read.table(parameter_df$cell_anno_file[i], header = T, sep="\t")
gene_anno <- read.table(parameter_df$gene_anno_file[i], header = T)
replace_mat <- read.table(parameter_df$perturbed_mat_file[i], header = T) %>% as.matrix()

#modify expression matrix
new_cell_anno <- rep("target_cell_type", 100) %>% set_names(cell_anno$cellid[cell_anno$is_target])
colData(data_sce)[["cell_ontology_class"]][match(names(new_cell_anno), colnames(data_sce))] = new_cell_anno


replace_mat <- out_mat_new
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

#test if gene annotation is right
