# ******************************************
# Analysis of the Tabula sapiens dataset
# ******************************************

if (!require("here")) {
  install.packages("here")
  library("here")
}

if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}

if (!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}

if (!require("scran")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("scran")
  library("scran")
}

if (!require("seismicGWAS")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

ts_obj = readRDS(here("all_data", "processed_expression", "TS_processed.rds"))


# the following lines have already been calculated for the processed data
# uncomment and run these when starting with raw data

# # quality control, normalization and label cleaning
# ts_obj$mt_counts <- assay(ts_obj, "counts")[grep(pattern = "^MT-", rowData(ts_obj)$gene_symbol), ] %>% colSums()
# ts_obj$mt_ratio <- ts_obj$mt_counts/ts_obj$n_counts_UMIs 
# ts_obj <- ts_obj[, ts_obj$mt_ratio<=0.1]

# # filter genes
# ts_obj <- ts_obj[(rowSums( assay(ts_obj, "counts")>0) >=10 &  rowSums( assay(ts_obj, "counts"))>=20),]

# # normalization
# cluster = quickCluster(ts_obj,  assay.type = "counts") 
# size.factor = calculateSumFactors(ts_obj, cluster=cluster, min.mean=0.1,assay.type="counts")
# ts_obj = logNormCounts(ts_obj, size.factors = size.factor,assay.type = "counts")


# calculate specificity score and cell type associations
ts_sscore  <- calc_specificity(sce = ts_obj , ct_label_col = "cluster_name", min_avg_exp_ct = 0.01)

# map to human genes
ts_sscore_hsa <- translate_gene_ids(ts_sscore, from = "hsa_ensembl")

# enrichment
magma_zscore_file <- list.files(here("all_data", "processed_gwas", "tm_gwas", "zscore"), full.names = T) %>%
  set_names(str_extract(string = ., pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.genes\\.out", replacement = "", x=.))

ts_association <- magma_zscore_file %>%
  map(~get_ct_trait_associations(sscore =ts_sscore_hsa, magma = .x))

ts_res <- ts_association %>%
  map(~select(.x, cell_type, pvalue)) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by="cell_type"))

write.table(ts_res, here("results" ,"Tabula_sapiens", "res.txt"),quote=F, sep="\t", row.names = F)
