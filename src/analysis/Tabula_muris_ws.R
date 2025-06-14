# ***************************************
# Tabula Muris window size analysis
#
# Note: MAGMA files of varying window
# sizes need to be first generated with
# tools/magma_gene_zscore_analysis.sh
# ***************************************
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

# load the object that's already processed in Tabula_muris_analysis.R
load(here("all_data", "processed_expression", "TM_processed.rda"))

# calculate specificity score and map to human genes
facs_sscore  <- calc_specificity(sce = facs_obj_sce , ct_label_col = "cluster_name")

# map to human genes
facs_sscore_hsa <- translate_gene_ids(facs_sscore, from = "mmu_symbol")

# varying window size
path_to_window_size_files 'data/gwas/zscore_varied_ws' # replace this with the path MAGMA run with various window sizes
gwas_file <- list.files(here(path_to_window_size_files)) %>%
  set_names(unlist(.)) %>% 
  map(~list.files(here(path_to_window_size_files,.x))) %>%
  map2(names(.), ~here(path_to_window_size_files,.y,.x)) %>%
  map(~set_names(.x, sub(pattern =".*\\/([^\\/]+)\\.genes\\.out$" , replacement = "\\1", x=.x)))

facs_res_ws <- gwas_file %>%
  map(~map(.x, ~get_ct_trait_associations(sscore = facs_sscore_hsa, magma = .x))) %>%
  map(~map(.x, ~select(.x, cell_type, pvalue))) %>%
  map(~map2(.x, names(.x), ~set_colnames(.x, c("cell_type", .y)))) %>%
  map(~purrr::reduce(.x, ~left_join(.x, .y, by="cell_type")))