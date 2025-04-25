library(SingleCellExperiment)
library(tidyverse)
library(here)
library(magrittr)

##### generate random seed for real data ####
load(here("data","expr","Tabula_muris","facs_clean.rda"))

#random sample cell types
facs_obj_sce <- facs_obj_sce[,!is.na(facs_obj_sce$cell_ontology_class)]
facs_obj_sce$cellid <- paste0("cell_", 1:ncol(facs_obj_sce))
colnames(facs_obj_sce) <- facs_obj_sce$cellid

set.seed(101)
random_ct <- sample(unique(facs_obj_sce$cell_ontology_class), 10)

##### sample cells according to ontology
random_ct_subset <- colData(facs_obj_sce) %>%
  as.data.frame() %>%
  add_count(cell_ontology_class, name = "cell_number") %>%
  distinct(cell_ontology_class, tissue, cell_number) %>%
  add_count(cell_ontology_class) %>%
  filter(n == 1) %>%
  filter(cell_ontology_class %in% random_ct)

#10 random cells
extra_rand_cell_type <- colData(facs_obj_sce) %>%
  as.data.frame() %>%
  add_count(cell_ontology_class, name = "cell_number") %>%
  distinct(cell_ontology_class, tissue, cell_number) %>%
  add_count(cell_ontology_class) %>%
  filter(n == 1) %>%
  filter(!cell_ontology_class %in% random_ct) %>%
  slice_sample(n=1)

#contactnate
random_ct_subset <- random_ct_subset %>% rbind(extra_rand_cell_type)

#create new parameter df 
cell_ratio <- c(0.05, 0.1, 0.25, 0.5, 0.75, 1.0) 

new_para_df <- expand.grid(cell_type = random_ct_subset$cell_ontology_class, cell_ratio = cell_ratio) %>%
  mutate(tissue = random_ct_subset$tissue[match(cell_type, random_ct_subset$cell_ontology_class)], 
         original_cell_number = random_ct_subset$cell_number[match(cell_type, random_ct_subset$cell_ontology_class)]) %>%
  mutate(cell_type = as.character(cell_type))

other_idx <- map2(new_para_df$cell_type, new_para_df$tissue,  ~facs_obj_sce$cellid[facs_obj_sce$cell_ontology_class != .x & facs_obj_sce$tissue == .y]) %>%
  map2(ceiling(new_para_df$cell_ratio * new_para_df$original_cell_number), ~ {all_idx = .x; num = .y; map(1:5, ~sample(all_idx, size = num))}) %>%
  map(~map2(.x, 1:length(.x), ~data.frame( .x) %>% set_colnames(paste0("sample_", .y)))) %>%
  map(~purrr::reduce(.x, ~cbind(.x, .y)))

#reset random cell types
random_ct[!random_ct %in% random_ct_subset$cell_ontology_class] = extra_rand_cell_type$cell_ontology_class

#write out cell index
new_para_df <- new_para_df %>% 
  mutate(cell_type_no = paste0("cell_type_", as.numeric(factor(cell_type, levels =random_ct)))) %>%
  mutate(cell_ratio = paste0("cr_", cell_ratio)) %>%
  mutate(output_header = here("results", "score_robustness", "cell_type_sample_by_tissue", cell_ratio, cell_type_no)) %>%
  mutate(seed_sample_file = here("data", "expr", "score_robustness", "cell_type_sample_by_tissue", paste0(cell_ratio, ".", cell_type_no, ".", "sample.txt")) )

new_para_df %>% write.table(here("data", "expr", "score_robustness", "cell_type_sample_by_tissue", "new_para_df.txt"), row.names = F, col.names = T, quote = F, sep="\t") 

pmap(list(other_idx, new_para_df$cell_ratio, new_para_df$cell_type_no), 
     ~write.table(..1, file = here("data", "expr", "score_robustness", "cell_type_sample_by_tissue", paste0(..2, ".", ..3, ".", "sample.txt")), row.names = F, col.names = T, quote = F,sep = "\t"))

