if (!require("here")) {
  install.packages("here")
  library("here")
}
if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}
if(!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("SingleCellExperiment")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SingleCellExperiment")
  library("SingleCellExperiment")
} 


##### load objects for standard fitted scDesign3 model ####
sc_para_file <- list.files(here("data","expr","causal_sim","standard_scdesign3_model"), pattern="scdesign3_para.rda", full.names = T) %>%
  set_names(str_extract(., pattern = "expr_ds_[0-9]*")) 

sc_para_list <- sc_para_file %>%
  map(~ {
    load(.x)
    list(copula = sce_copula,
         para = sce_para)
  }) 

save(sc_para_list, file = here("data","expr","causal_sim","standard_scdesign3_model","sc_para_list.rda"))


##### load objects for fitted scDeisn3 model with 3 cell types ####
sc_para_file <- list.files(here("data","expr","causal_sim","multi_ct_scdesign3_model"), pattern="scdesign3_para.rda", full.names = T) %>%
  set_names(str_extract(., pattern = "expr_ds_[0-9]*")) 

sc_para_list <- sc_para_file %>%
  map(~ {
    load(.x)
    list(copula = sce_copula,
         para = sce_para)
  }) 

save(sc_para_list, file = here("data","expr","causal_sim","multi_ct_scdesign3_model","scdesign3_para.rda"))
