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

#redo analysis
data("mmu_hsa_mapping") 
facs_obj_sce = reset_seismic_analysis(facs_obj_sce)
facs_obj_sce  = cal_stat(data_obj = facs_obj_sce , meta_data = as.data.frame(colData(facs_obj_sce)), group = "cluster_name")
facs_obj_sce  = cal_sscore(data_obj = facs_obj_sce) 
facs_obj_sce = trans_mmu_to_hsa_stat(facs_obj_sce , gene_mapping_table=mmu_hsa_mapping, from="mmu_symbol", to="hsa_entrez")
facs_obj_sce = add_glob_stats(facs_obj_sce, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 

#varying window size
gwas_zscore = list.files(here("data","gwas","tm_gwas","zscore_r")) %>%
  set_names(unlist(.)) %>% 
  map(~load_zscore(here("data","gwas","tm_gwas","zscore_varied_ws",.x))) 

facs_res_ws = gwas_zscore %>%
  map(~cal_ct_asso(facs_obj_sce, .x, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.1& max_exp_ct>0.1", asso_model = "linear")) %>%
  map(~get_ct_asso(.x, trait_name = "all", asso_model = "linear", merge_output = T))

#save results
map2(facs_res_ws, names(facs_res_ws), ~write.table(.x, here("results","varied_ws",paste0("all_res.",.y,".txt")),quote=F, sep="\t", row.names = F))
