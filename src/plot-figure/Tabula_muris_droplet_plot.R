#generate plots for results of the Tabula Muris dataset
##### load packages and results from scdrs/ours/fuma/magma #####
if (!require("here")){
  install.packages("here")
  library("here")
}
if (!require("tidyverse")){
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("magrittr")){
  install.packages("magrittr")
  library("magrittr")
}
if (!require("here")){
  install.packages("here")
  library("here")
}

#color mappings
color_mapping_vec = c("scDRS" ="#46B8DAFF" , "seismic" = "#D43F3AFF" , "S-MAGMA"="#EEA236FF" , "FUMA"="#5CB85CFF")

#load results - droplet
droplet_res = read.table( here("results","Tabula_muris","droplet","seismic","droplet_res.txt"), header = T, sep = "\t") %>% as_tibble()

droplet_ct_meta = tibble(cell_type = as.character(droplet_res$cell_type)) %>%
  mutate(tissue = strsplit(cell_type, split=".", fixed = T) %>% map(~.x[1]) %>% unlist) %>%
  mutate(cell_ontology = strsplit(cell_type, split=".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(tissue = gsub(pattern = "Heart_and_|Limb_|_Gland",replacement = "", x = tissue) ) %>%
  mutate(cell_ontology = gsub(pattern = "alveolar epithelial type 1 cells, alveolar epithelial type 2 cells, club cells, and basal cells", "alveolar epithelial cells",x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern = "dendritic cells, alveolar macrophages, and interstital macrophages","dendritic cells and alveolar macrophages",x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern=" positive","+",x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern = " negative", "-", x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern = "Thymocytes", "T cells", x=cell_ontology,ignore.case=T)) %>%
  mutate(tissue_new = tissue) %>%
  mutate(tissue_new = ifelse(tissue_new == "Marrow","Blood/Immune",tissue_new)) %>%
  mutate(tissue_new = ifelse(grepl(pattern = "B cell|T cell|immune cells|erythrocyte|macrophage|microglia|myeloid|leukocyte|natural killer|antigen|NK|monocyte|lymphocyte|blood|Macrophage|thymocyte", x=cell_ontology), "Blood/Immune",tissue_new)) %>%
  mutate(ontology_new = ifelse(tissue_new == "Blood/Immune", paste0(tissue,".",cell_ontology),cell_ontology)) 
