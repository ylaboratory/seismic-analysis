#generate plots for results of the Tabula sapiens data set
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

#load results - ts
ts_res <- read.table( here("results","Tabula_sapiens","seismic","new_ts_res.txt"), header = T, sep = "\t") %>%
  as_tibble() %>%
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))  %>%
  mutate(cell_type = factor(cell_type, levels = sort(cell_type))) %>%
  arrange(cell_type)

#load results - other frameworks
ts_scdrs_res <- list.files(here("results","Tabula_sapiens","scDRS"),pattern = "*.scdrs_group.cluster_name",full.names = T) %>% 
  set_names(str_extract(string = ., pattern = "(?<=/)[^/]+$") %>% gsub(pattern = ".scdrs_group.cluster_name", replacement = "", x=.)) %>%
  map(~read.table(.x, header=T, sep="\t") %>% as_tibble()) %>%
  map(~mutate(.x, P_val = pnorm(assoc_mcz, lower.tail = F))) %>%
  map(~select(.x, any_of(c("group","P_val")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>% 
  purrr::reduce(~left_join(.x, .y, by = "cell_type")) %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))%>%
  mutate(cell_type = factor(cell_type, levels = levels(ts_res$cell_type))) %>%
  arrange(cell_type)

ts_fuma_res <- list.files(here("results","new_Tabula_sapiens","FUMA"), pattern = "*.gsa.out", full.names = T) %>%
  set_names(str_extract(string = . , pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.[A-Za-z.]+$", replacement = "", x= .)) %>%
  map(~read.table(.x, header = T) %>% as_tibble()) %>%
  map(~select(.x, any_of(c("VARIABLE", "P")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by = "cell_type"))  %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) 

ts_magma_res <- list.files(here("results","new_Tabula_sapiens","S-MAGMA"), pattern = "*.gsa.out", full.names = T) %>%
  set_names(str_extract(string = . , pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.[A-Za-z.]+$", replacement = "", x= .)) %>%
  map(~read.table(.x, header = T) %>% as_tibble()) %>%
  map(~select(.x, any_of(c("VARIABLE", "P")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by = "cell_type"))  %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) 

#encode cell types
ts_ct_meta <- tibble(cell_type = as.character(ts_res$cell_type)) %>%
  mutate(tissue = strsplit(cell_type, split=".", fixed = T) %>% map(~.x[1]) %>% unlist) %>%
  mutate(cell_ontology = strsplit(cell_type, split=".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(cell_ontology = gsub(pattern=" positive|-positive","+",x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern = " negative|-negative", "-", x=cell_ontology)) %>% 
  mutate(tissue_new = ifelse(tissue=="Bone_Marrow","Blood/Immune",tissue)) %>%
  mutate(tissue_new = ifelse(grepl(pattern = "B cell|T cell|immune cells|erythrocyte|macrophage|microglia|myeloid|leukocyte|natural killer|antigen|NK|monocyte|lymphocyte|blood|thymocyte", x=cell_ontology), "Blood/Immune",tissue_new)) %>%
  mutate(ontology_new = ifelse(tissue_new == "Blood/Immune", paste0(tissue,".",cell_ontology),cell_ontology)) 

#harmonise cell types for FUMA/S-MAGMA output
ts_fuma_ct_mapping <- read.table(here("data","expr","Tabula_sapiens","new_ts.magma.aux.txt"), header = T, sep = "\t") %>% as_tibble()
ts_magma_ct_mapping <- read.table(here("data","expr","Tabula_sapiens","new_ts.fuma.aux.txt"), header = T, sep = "\t") %>% as_tibble()

ts_fuma_res <- ts_fuma_ct_mapping %>% 
  left_join(ts_fuma_res, by = c("encoded_name" = "cell_type")) %>%
  select( -encoded_name )  %>% 
  mutate(cell_type = factor(cell_type, levels = levels(ts_res$cell_type))) %>%
  arrange(cell_type)

ts_magma_res <- ts_magma_ct_mapping %>% 
  left_join(ts_magma_res, by = c("encoded_name" = "cell_type")) %>%
  select( -encoded_name ) %>% 
  mutate(cell_type = factor(cell_type, levels = levels(ts_res$cell_type))) %>%
  arrange(cell_type)

#### heatmap ####
library(ComplexHeatmap)
library(circlize)
ht_opt(TITLE_PADDING=unit(5,"mm"))
col_fun <- colorRamp2(c(-0.5, -log10(0.05), 4,10), viridis::viridis(4))

ts_fdr_mat <- ts_res %>%
  mutate(across(where(is.double), ~p.adjust(.x, method ="fdr"))) %>%
  select(all_of(c("cell_type",neuropsy_disesaes,immune_diseases,others))) %>%
  select(where(is.double)) %>%
  as.matrix() %>%
  t() %>%
  set_colnames( ts_ct_meta$cell_type[match(ts_res$cell_type, ts_ct_meta$cell_type)]) %>%
  set_rownames(trait_meta$official_names[match(rownames(.), trait_meta$trait_names)])

Heatmap(-log10(ts_fdr_mat), name = "-log10(FDR)",cluster_rows = F,cluster_columns = F,show_column_dend = F, show_row_dend = F,rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
        column_split = ts_ct_meta$tissue_new[match(colnames(ts_fdr_mat), ts_ct_meta$cell_type)],
        column_labels = ts_ct_meta$ontology_new[match(colnames(ts_fdr_mat), ts_ct_meta$cell_type)],
        row_split = trait_meta$trait_type[match(rownames(ts_fdr_mat), trait_meta$official_names)] %>% factor(levels = c("neuropsy","immune","others")),
        row_names_side = "left",
        column_title_gp = gpar(fontsize=11),
        row_names_gp = gpar(fontsize=9),
        column_names_gp = gpar(fontsize=8),
        column_names_rot = 55,
        column_title_rot = 30,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(ts_fdr_mat[i, j] < 1e-2) {
            grid.text("**", x, y, gp = gpar(fontsize = 8))
          } else if(ts_fdr_mat[i, j] <= 0.05) {
            grid.text("*", x, y, gp = gpar(fontsize = 8))
          }
        })

#### correlation ####
all_results_ts <- list(seismic = ts_res, scDRS = ts_scdrs_res, FUMA = ts_fuma_res, `S-MAGMA` = ts_magma_res) %>%
  map(~pivot_longer(.x, !cell_type, names_to = "trait",values_to = "Pvalue")) %>% 
  map2(names(.), ~mutate(.x, method = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>%
  group_by(trait, method) %>% 
  mutate(fdr = p.adjust(Pvalue, method="fdr")) %>%
  ungroup() %>%
  mutate(trait_name = trait_meta$official_names[match(trait,trait_meta$trait_names)]) %>%
  mutate(trait_type = trait_meta$trait_type[match(trait,trait_meta$trait_names)]) %>%
  select(cell_type,trait_name, trait_type, method, Pvalue, fdr)

cor_ts_metric <- all_results_ts %>%
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
  select(trait_name, between, correlation) %>%
  mutate(is_negative = ifelse(correlation<0, "negative", "positive"))

ggplot( cor_ts_metric, aes(x=between, y=correlation, color=is_negative )) + 
  geom_point(size=2.5) +
  geom_segment(aes(xend=between, yend=0), linetype="solid",linewidth=0.75) +
  scale_color_manual(name="value",values =c("grey","#2B3990FF") ) +
  theme_classic() + 
  ylim(-0.3,1) + 
  coord_flip() + 
  facet_wrap(~trait_name,nrow=5) + 
  ggtitle("TM ts data set, Spearman correlation across all cell types") +
  ylab("Spearman's correlation") + 
  theme(plot.title = element_text(hjust = 0.5))
