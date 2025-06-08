#generate plots for results of the Tabula Muris droplet dataset
##### load packages and results from scdrs/ours/fuma/magma #####
if (!require("here")) {
  install.packages("here")
  library("here")
}
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}
if (!require("here")) {
  install.packages("here")
  library("here")
}

#color mappings
color_mapping_vec <- c("scDRS" ="#46B8DAFF" , "seismic" = "#D43F3AFF" , "S-MAGMA"="#EEA236FF" , "FUMA"="#5CB85CFF")

#encode trait mapping
neuropsy_disesaes <- c("insominia","BD","MDD","Nrt_new","Scz_new")
immune_diseases <- c("AID","Hypothyroidism","RA","ibd","cd","uc")
others <- c("Smoking","BMI","College_edu","RBC","Lymphocyte_count","Monocyte_count","T1D","T2D_2","Cardiovas","SBP","AF","CKD","glucose_q","HDL_q","LDL_q","TG_q")
trait_meta <- tibble(trait_names = c(neuropsy_disesaes, immune_diseases, others),
                    official_names = c("Insomnia","Bipolar disorder","Depression","Neuroticism","Schizophrenia",
                                       "Autoimmune diseases","Hypothroidism","Rheumatoid arthritis","Inflammatory bowel disease","Crohn's disease","Ulcerative colitis",
                                       "Smoking","BMI","College education","Erythrocyte count","Lymphocyte count","Monocyte count","Type I diabetes","Type II diabetes","Cardiovascular diseases",
                                       "Systolic blood pressure","Atrial fibrillation","Chronic kidney disease","Glucose level","HDL level","LDL level","Triglyceride level"),
                    trait_type = factor(c(rep("neuropsy",length(neuropsy_disesaes)), rep("immune",length(immune_diseases)), rep("others",length(others))), levels=c("neuropsy","immune","others")) )

#load results - droplet
droplet_res <- read.table( here("results","Tabula_muris","droplet","seismic","new_droplet_res.txt"), header = T, sep = "\t") %>%
  as_tibble() %>%
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))  %>%
  mutate(cell_type = factor(cell_type, levels = sort(cell_type))) %>%
  arrange(cell_type)

#load results - other frameworks
droplet_scdrs_res <- list.files(here("results","Tabula_muris","droplet","scDRS"),pattern = "*.scdrs_group.cluster_name",full.names = T) %>% 
  set_names(str_extract(string = ., pattern = "(?<=/)[^/]+$") %>% gsub(pattern = ".scdrs_group.cluster_name", replacement = "", x=.)) %>%
  map(~read.table(.x, header=T, sep="\t") %>% as_tibble()) %>%
  map(~mutate(.x, P_val = pnorm(assoc_mcz, lower.tail = F))) %>%
  map(~select(.x, any_of(c("group","P_val")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>% 
  purrr::reduce(~left_join(.x, .y, by = "cell_type")) %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))%>%
  mutate(cell_type = factor(cell_type, levels = levels(droplet_res$cell_type))) %>%
  arrange(cell_type)

droplet_fuma_res <- list.files(here("results","Tabula_muris","droplet","FUMA"), pattern = "*.gsa.out", full.names = T) %>%
  set_names(str_extract(string = . , pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.[A-Za-z.]+$", replacement = "", x= .)) %>%
  map(~read.table(.x, header = T) %>% as_tibble()) %>%
  map(~select(.x, any_of(c("VARIABLE", "P")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by = "cell_type"))  %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) 

droplet_magma_res <- list.files(here("results","Tabula_muris","droplet","S-MAGMA"), pattern = "*.gsa.out", full.names = T) %>%
  set_names(str_extract(string = . , pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.[A-Za-z.]+$", replacement = "", x= .)) %>%
  map(~read.table(.x, header = T) %>% as_tibble()) %>%
  map(~select(.x, any_of(c("VARIABLE", "P")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by = "cell_type"))  %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) 

droplet_ct_meta <- tibble(cell_type = as.character(droplet_res$cell_type)) %>%
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

#harmonise cell types for FUMA/S-MAGMA output
droplet_fuma_ct_mapping <- read.table(here("data","expr","Tabula_muris","new_tm_droplet.magma.aux.txt"), header = T, sep = "\t") %>% as_tibble()
droplet_magma_ct_mapping <- read.table(here("data","expr","Tabula_muris","new_tm_droplet.fuma.aux.txt"), header = T, sep = "\t") %>% as_tibble()

droplet_fuma_res <- droplet_fuma_ct_mapping %>% 
  left_join(droplet_fuma_res, by = c("encoded_name" = "cell_type")) %>%
  select( -encoded_name )  %>% 
  mutate(cell_type = factor(cell_type, levels = levels(droplet_res$cell_type))) %>%
  arrange(cell_type)

droplet_magma_res <- droplet_magma_ct_mapping %>% 
  left_join(droplet_magma_res, by = c("encoded_name" = "cell_type")) %>%
  select( -encoded_name ) %>% 
  mutate(cell_type = factor(cell_type, levels = levels(droplet_res$cell_type))) %>%
  arrange(cell_type)

#### heatmap ####
library(ComplexHeatmap)
library(circlize)
ht_opt(TITLE_PADDING=unit(5,"mm"))
col_fun <- colorRamp2(c(-0.5, -log10(0.05), 4,10), viridis::viridis(4))

droplet_fdr_mat <- droplet_res %>%
  mutate(across(where(is.double), ~p.adjust(.x, method ="fdr"))) %>%
  select(all_of(c("cell_type",neuropsy_disesaes,immune_diseases,others))) %>%
  select(where(is.double)) %>%
  as.matrix() %>%
  t() %>%
  set_colnames( droplet_ct_meta$cell_type[match(droplet_res$cell_type, droplet_ct_meta$cell_type)]) %>%
  set_rownames(trait_meta$official_names[match(rownames(.), trait_meta$trait_names)])

Heatmap(-log10(droplet_fdr_mat), name = "-log10(FDR)",cluster_rows = F,cluster_columns = F,show_column_dend = F, show_row_dend = F,rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
        column_split = droplet_ct_meta$tissue_new[match(colnames(droplet_fdr_mat), droplet_ct_meta$cell_type)],
        column_labels = droplet_ct_meta$ontology_new[match(colnames(droplet_fdr_mat), droplet_ct_meta$cell_type)],
        row_split = trait_meta$trait_type[match(rownames(droplet_fdr_mat), trait_meta$official_names)] %>% factor(levels = c("neuropsy","immune","others")),
        row_names_side = "left",
        column_title_gp = gpar(fontsize=11),
        row_names_gp = gpar(fontsize=9),
        column_names_gp = gpar(fontsize=8),
        column_names_rot = 55,
        column_title_rot = 30,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(droplet_fdr_mat[i, j] < 1e-2) {
            grid.text("**", x, y, gp = gpar(fontsize = 8))
          } else if(droplet_fdr_mat[i, j] <= 0.05) {
            grid.text("*", x, y, gp = gpar(fontsize = 8))
          }
        })

#### correlation ####
all_results_droplet <- list(seismic = droplet_res, scDRS = droplet_scdrs_res, FUMA = droplet_fuma_res, `S-MAGMA` = droplet_magma_res) %>%
  map(~pivot_longer(.x, !cell_type, names_to = "trait",values_to = "Pvalue")) %>% 
  map2(names(.), ~mutate(.x, method = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>%
  group_by(trait, method) %>% 
  mutate(fdr = p.adjust(Pvalue, method="fdr")) %>%
  ungroup() %>%
  mutate(trait_name = trait_meta$official_names[match(trait,trait_meta$trait_names)]) %>%
  mutate(trait_type = trait_meta$trait_type[match(trait,trait_meta$trait_names)]) %>%
  select(cell_type,trait_name, trait_type, method, Pvalue, fdr)

cor_droplet_metric <- all_results_droplet %>%
  group_by(trait_name, method) %>%
  mutate(method = factor(method, levels=c("S-MAGMA","FUMA","scDRS","seismic"))) %>% 
  arrange(cell_type) %>%
  summarize(all_p = list(Pvalue)) %>%
  inner_join(.,., by="trait_name",relationship = "many-to-many") %>%
  filter(as.numeric(method.x) > as.numeric(method.y)) %>%
  rowwise() %>%
  mutate(correlation = cor(unlist(all_p.x), unlist(all_p.y), method="spearman")) %>%
  mutate(between = paste0(method.x, " and ",method.y)) %>%
  mutate(between = factor(between, levels =rev(c("seismic and scDRS","seismic and FUMA","seismic and S-MAGMA","scDRS and FUMA","scDRS and S-MAGMA","FUMA and S-MAGMA")))) %>%
  select(trait_name, between, correlation) 

ggplot(cor_droplet_metric, aes(x=between, y=correlation)) + 
  geom_point(size=2.5, color="#2B3990FF") +
  geom_segment(aes(xend=between, yend=0), linetype="solid",linewidth=0.75, color="#2B3990FF") +
  theme_classic() + 
  ylim(-0.1,1) + 
  coord_flip() + 
  facet_wrap(~trait_name,nrow=5) + 
  ggtitle("TM droplet data set, Spearman correlation across all cell types") +
  ylab("Spearman's correlation") + 
  theme(plot.title = element_text(hjust = 0.5))
