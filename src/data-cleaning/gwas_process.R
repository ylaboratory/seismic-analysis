#script to generate data for causal simulation
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

#load objects
load(here("data","ref","mapping","mmu_hsa_mapping.rda"))

z_score_file = list.files(here("data","gwas","tm_gwas","zscore"), full.names = T) %>% 
  set_names(str_extract(. , pattern = "(?<=/)[^/]+$" ) %>% gsub(pattern = ".genes.out",replacement = "", fixed = T )) %>%
  .[match(unique(names(.)), names(.))]
z_stat = z_score_file %>% 
  map(~read.table(.x, header=T)%>% as_tibble()) 

z_stat = z_stat %>% 
  map(~left_join(.x, mmu_hsa_mapping %>% distinct(hsa_symbol,hsa_entrez), by = c("GENE"="hsa_entrez"))) %>%
  map(~distinct(.x)) %>%
  map(~drop_na(.x,hsa_symbol)) %>%
  map(~select(.x, hsa_symbol, ZSTAT)) %>%
  map(~group_by(.x, hsa_symbol)) %>%
  map(~mutate(.x,ZSTAT = mean(ZSTAT))) %>%
  map(~ungroup(.x)) %>%
  map(~distinct(.x,  hsa_symbol,ZSTAT)) %>% 
  map2(names(z_score_file), ~set_colnames(.x, c("GENE",.y)))

z_stat %>% map2(names(z_stat), ~write_tsv(.x, file=paste0("data/gwas/tm_gwas/scdrs_gs/",.y,".tsv")))

#for neuron gs
z_score_file_neuron = list.files(here("data","gwas","brain_gwas","zscore"), full.names = T) %>% 
  set_names(str_extract(. , pattern = "(?<=/)[^/]+$" ) %>% gsub(pattern = ".genes.out",replacement = "", fixed = T )) %>%
  .[match(unique(names(.)), names(.))]
z_stat_neuron = z_score_file_neuron %>% 
  map(~read.table(.x, header=T)%>% as_tibble()) 

z_stat_neuron = z_stat_neuron %>% 
  map(~left_join(.x, mmu_hsa_mapping %>% distinct(hsa_symbol,hsa_entrez), by = c("GENE"="hsa_entrez"))) %>%
  map(~distinct(.x)) %>%
  map(~drop_na(.x,hsa_symbol)) %>%
  map(~select(.x, hsa_symbol, ZSTAT)) %>%
  map(~group_by(.x, hsa_symbol)) %>%
  map(~mutate(.x,ZSTAT = mean(ZSTAT))) %>%
  map(~ungroup(.x)) %>%
  map(~distinct(.x,  hsa_symbol,ZSTAT)) %>% 
  map2(names(z_score_file_neuron), ~set_colnames(.x, c("GENE",.y)))

z_stat_neuron %>% map2(names(z_stat_neuron), ~write_tsv(.x, file=paste0("data/gwas/brain_gwas/scdrs_gs/",.y,".tsv")))
