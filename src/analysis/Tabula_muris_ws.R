#Analysis of Tabula muris data set
##### 1. load packages and data#######
###load packages
if (!require("here")){
  install.packages("here")
  library("here")
}

if (!require("magrittr")){
  install.packages("magrittr")
  library("magrittr")
}
if (!require("tidyverse")){
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("scran")){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("scran")
  library("scran")
}  #normalize data

if (!require("seismicGWAS")){
  if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

#load the object that's filtered already in the Tabula_muris_analysis scripts
load(here("data","expr","Tabula_muris","TM_processed.rda"))

#calculate specificity score and map to human genes
facs_sscore  <- calc_specificity(sce = facs_obj_sce , ct_label_col = "cluster_name")

#map to human genes
facs_sscore_hsa <- translate_gene_ids(facs_sscore, from = "mmu_symbol")

#varying window size
gwas_file <- list.files(here("data","gwas","tm_gwas","zscore_varied_ws")) %>%
  set_names(unlist(.)) %>% 
  map(~list.files(here("data","gwas","tm_gwas","zscore_varied_ws",.x))) %>%
  map2(names(.), ~here("data","gwas","tm_gwas","zscore_varied_ws",.y,.x)) %>%
  map(~set_names(.x, sub(pattern =".*\\/([^\\/]+)\\.genes\\.out$" , replacement = "\\1", x=.x)))

facs_res_ws <- gwas_file %>%
  map(~map(.x, ~get_ct_trait_associations(sscore = facs_sscore_hsa, magma = .x))) %>%
  map(~map(.x, ~select(.x, cell_type, pvalue))) %>%
  map(~map2(.x, names(.x), ~set_colnames(.x, c("cell_type", .y)))) %>%
  map(~purrr::reduce(.x, ~left_join(.x, .y, by="cell_type")))

#save results
map2(facs_res_ws, names(facs_res_ws), ~write.table(.x, here("results","varied_ws",paste0("new_all_res.",.y,".txt")),quote=F, sep="\t", row.names = F))
