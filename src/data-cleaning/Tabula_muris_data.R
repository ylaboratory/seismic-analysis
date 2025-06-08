#Preprocessing steps for Tabula muris data set
##### 1. load packages and data#######
### 1 load packages
if(!require("Seurat")) {
  install.packages("Seurat")
  library("Seurat")
}
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
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("SingleCellExperiment")
  library("SingleCellExperiment")
}  #normalize data

###2 load data
facs_obj_file = list.files(path = here("raw","expr","Tabula_muris","Seurat"), pattern = "facs.*Robj")
droplet_obj_file = list.files(path = here("raw","expr","Tabula_muris","Seurat"), pattern = "droplet.*Robj")

##### 2. merge data#######
facs_obj = facs_obj_file %>%
  map(~{load(here("raw","expr","Tabula_muris","Seurat",.x)); tiss}) %>%
  map(~UpdateSeuratObject(.x)) %>%
  reduce(~merge(.x,.y))
droplet_obj = droplet_obj_file %>%
  map(~{load(here("raw","expr","Tabula_muris","Seurat",.x)); tiss}) %>%
  map(~UpdateSeuratObject(.x)) %>%
  reduce(~merge(.x,.y))

##### 3. transform data to sce and clean labels #####
facs_obj_sce = as.SingleCellExperiment(facs_obj)
droplet_obj_sce = as.SingleCellExperiment(droplet_obj)

#clean cell ontology label and remove cells without such labels
colData(facs_obj_sce) = colData(facs_obj_sce) %>% 
  as_tibble() %>% #easy to process
  group_by(cell_ontology_id) %>% #unify cell_ontology_class name (some typo)
  mutate(cell_ontology_class = cell_ontology_class[1]) %>%
  ungroup %>%
  DataFrame()

colData(droplet_obj_sce) = colData(droplet_obj_sce) %>% 
  as_tibble() %>% #easy to process
  group_by(cell_ontology_id) %>% #unify cell_ontology_class name (some typo)
  mutate(free_annotation = ifelse(free_annotation == "DN4-DP in transition Cd69_negative thymocytes","DN4-DP in transition Cd69 negative thymocytes" ,free_annotation)) %>%
  mutate(free_annotation = ifelse(free_annotation == "DN4-DP transition Cd69 negative rapidly dividing thymocytes","DN4-DP in transition Cd69 negative rapidly dividing thymocytes" ,free_annotation)) %>%
  mutate(cell_ontology_class = cell_ontology_class[1]) %>%
  ungroup %>%
  DataFrame()

save(facs_obj_sce, file = here("data","expr","Tabula_muris","facs_clean.rda"))
save(droplet_obj_sce, file = here("data","expr","Tabula_muris","droplet_clean.rda"))