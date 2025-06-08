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
if (!require("seismicGWAS")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}
library(ComplexHeatmap)
library(circlize)
library(cowplot)

source(here("src","tools","magma_fuma_file_prep.R"))
source(here("src","tools","sparse_mat_util.R"))

#### calculate specificity score of Tabula muris datasets ####
load(here("data","expr","Tabula_muris","TM_processed.rda"))

###specificity score 
facs_sscore  <- calc_specificity(sce = facs_obj_sce , ct_label_col = "cluster_name")
droplet_sscore  <- calc_specificity(sce = droplet_obj_sce , ct_label_col = "cluster_name", min_avg_exp_ct = 0.01) #for droplet


#### plot distribution for marker genes ####
#markers
pancreas_marker <- list("Pancreas.pancreatic A cell" = c("Irx2", "Ttr", "Gcg", "Slc7a2", "Nxph1"),
                          "Pancreas.beta cell" = c("Ins1", "Ins2", "G6pc2", "Iapp", "Spock2"),
                          "Pancreas.pancreatic D cell" = c( "Sst", "Spock3"),
                          "Pancreas.pancreatic PP cell" = c("Ppy"))


#plot marker distribution #1 heatmap 
submat_sscore_pancreas <- facs_sscore[match(unlist(pancreas_marker), rownames(facs_sscore)), match(names(pancreas_marker), colnames(facs_sscore))] %>%
  t()

#heatmap
col_fun <- colorRamp2(c(0,0.15), c("white", "navy"))

Heatmap(submat_sscore_pancreas, name = "seismic specificity score",cluster_rows = F,cluster_columns = F,show_column_dend = F, show_row_dend = F,
                        rect_gp = gpar(col = "grey", lwd = 0.5), 
                        col = col_fun, 
                        column_split = c(rep("Alpha cell markers", 5), rep("Beta cell markers", 5), rep("Delta cell markers", 2), "PP cell markers"),
                        row_names_side = "right",
                        column_names_rot = 55,
                        column_title_rot = 30,
                        row_title_rot = 60,
                        column_title_gp = gpar(fontsize=9),
                        row_title_gp = gpar(fontsize = 9),
                        row_names_gp = gpar(fontsize=8),
                        column_names_gp = gpar(fontsize=8))

#plot marker distribution
pancreas_sscore_rank <- tibble(gene = unlist(pancreas_marker),
                               target_cell_type = map2(names(pancreas_marker), map(pancreas_marker, ~length(.x)), ~rep(.x, .y)) %>% unlist)

pancreas_sscore_rank_all <- facs_sscore[match(unlist(pancreas_marker), rownames(facs_sscore)),] %>%
  t() %>%
  as_tibble(rownames = "cell_type") %>%
  pivot_longer(cols = !cell_type, names_to = "gene", values_to = "sscore") %>%
  group_by(gene) %>%
  arrange(-sscore) %>%
  filter(1:n() <= 20) %>%
  ungroup() %>%
  left_join(pancreas_sscore_rank, by = "gene") %>%
  split(.$target_cell_type) %>%
  map(~split(.x, .x$gene)) %>%
  map(~map(.x, ~mutate(.x, cell_type_label = ifelse(cell_type %in% names(pancreas_marker), cell_type, "others")))) %>%
  map(~map(.x, ~mutate(.x, cell_type = factor(cell_type, levels = cell_type)))) %>%
  map(~map(.x, ~mutate(.x, cell_type_label = factor(cell_type_label, levels = c(names(pancreas_marker), "others"))))) 

pancreas_sscore_rank_plot_list <- pancreas_sscore_rank_all %>%
  map(~map(.x, ~ggplot(.x, aes(x = cell_type, y = sscore, fill = cell_type_label)) + 
             geom_bar(stat = "identity", alpha = 0.85) + theme_classic() + 
             theme(axis.text.x = element_blank()) + 
             scale_fill_manual(values = ggsci::pal_npg()(5) %>% set_names(c(names(pancreas_marker), "others"))) +
             ylab(paste0("seismic specificity score")))) %>%
  map(~map2(.x,names(.x), ~.x + ggtitle(.y) + theme(plot.title = element_text(hjust = 0.5))))

empty_plot <- ggplot() + theme_void()

pancreas_sscore_rank_plot_list <- pancreas_sscore_rank_plot_list %>%
  map(~ {if (length(.x) == 5) .x else c(.x, list(empty_plot)[rep(1, 5 - length(.x))])}) %>%
  purrr::reduce(~c(.x, .y))

final_grid <- final_grid <- cowplot::plot_grid(plotlist = pancreas_sscore_rank_plot_list,ncol = 4, nrow = 5,align = "hv",byrow = F)

#plot violin plot comparison
hk_gene_list <- c("Actb", "Atp5f1", "B2m", "Gapdh", "Hprt", "Pgk1", "Rer1", "Rpl13a", "Rpl27", "Sdha", "Tbp", "Ubc")

vn_marker_df <- pancreas_sscore_rank_all %>% 
  map(~map(.x, ~filter(.x,cell_type == target_cell_type))) %>%
  map(~purrr::reduce(.x, ~rbind(.x, .y))) %>%
  purrr::reduce(~rbind(.x, .y)) %>% 
  select(-target_cell_type, -cell_type_label) %>%
  mutate(pair = "pancreatic cell type marker in target cell types") %>%
  rbind(facs_sscore[match(hk_gene_list, rownames(facs_sscore)),] %>% as_tibble(rownames = "cell_type") %>% 
          pivot_longer(cols = -cell_type, names_to = "gene", values_to = "sscore") %>%
          mutate(pair = "housekeeping genes in all cell types")) %>%
  rbind(facs_sscore[match(unlist(pancreas_marker), rownames(facs_sscore)), !colnames(facs_sscore) %in% names(pancreas_marker)] %>% 
          as_tibble(rownames = "cell_type") %>% 
          pivot_longer(cols = -cell_type, names_to = "gene", values_to = "sscore") %>%
          mutate(pair = "pancreatic cell type marker in non-pancreatic cell types"))

ggplot(vn_marker_df, aes(x = pair, y = sscore, fill = pair)) +
  # Boxplot with reduced outlier visibility
  geom_boxplot(
    width = 0.6, 
    alpha = 0.8,
    outlier.size = 0.7,     
    outlier.alpha = 0.3,    
    outlier.color = "grey40" 
  ) +
  ggsci::scale_fill_npg() + 
  labs(y = "seismic specificity score", x = "",fill = "Category") +
  theme_classic() +
  theme(
    legend.position = "right", legend.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_blank()
  )


#### heatmap between mean and sscore ####
facs_mean <- calc_ct_mean(facs_obj_sce, ct_label_col = "cluster_name")
droplet_mean <- calc_ct_mean(droplet_obj_sce, ct_label_col = "cluster_name")

#calculate correlation
facs_mean_cor <-  cor(as.matrix(t(facs_mean)))
facs_sscore_cor <- cor(as.matrix(facs_sscore))
droplet_mean_cor <-  cor(as.matrix(t(droplet_mean)))
droplet_sscore_cor <- cor(as.matrix(droplet_sscore))

#plot heatmap
#load cell type annotation
facs_ct_meta <- tibble(cell_type = colnames(facs_sscore)) %>%
  mutate(tissue = strsplit(cell_type, split=".", fixed = T) %>% map(~.x[1]) %>% unlist) %>%
  mutate(cell_ontology = strsplit(cell_type, split=".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(tissue = gsub(pattern = "_Non-Myeloid|_Myeloid|Large_|Limb_|_Gland",replacement = "", x = tissue) ) %>%
  mutate(cell_ontology = gsub(pattern = "alveolar epithelial type 1 cells, alveolar epithelial type 2 cells, club cells, and basal cells", "alveolar epithelial cells",x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern = "dendritic cells, alveolar macrophages, and interstital macrophages","dendritic cells and alveolar macrophages",x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern=" positive","+",x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern = " negative", "-", x=cell_ontology)) %>%
  mutate(cell_ontology = gsub(pattern = "Thymocytes", "T cells", x=cell_ontology,ignore.case=T)) %>% 
  mutate(cell_ontology = gsub(pattern = "thyomcyte", replacement = "thymocyte",x=cell_ontology)) %>%  #harmonised ontology names
  mutate(tissue_new = tissue) %>% 
  mutate(tissue_new = ifelse(tissue_new == "Marrow","Blood/Immune",tissue_new)) %>%
  mutate(tissue_new = ifelse(grepl(pattern = "B cell|T cell|immune cells|erythrocyte|macrophage|microglia|myeloid|leukocyte|natural killer|antigen|NK|monocyte|lymphocyte|blood|thymocyte|dendritic", x=cell_ontology), "Blood/Immune",tissue_new)) %>%
  mutate(ontology_new = ifelse(tissue_new == "Blood/Immune", paste0(tissue,".",cell_ontology),cell_ontology)) 

droplet_ct_meta <- tibble(cell_type = colnames(droplet_sscore)) %>%
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
  mutate(tissue_new = ifelse(grepl(pattern = "B cell|T cell|immune cells|erythrocyte|macrophage|microglia|myeloid|leukocyte|natural killer|antigen|NK|monocyte|lymphocyte|blood|Macrophage|thymocyte|dendritic", x=cell_ontology), "Blood/Immune",tissue_new)) %>%
  mutate(ontology_new = ifelse(tissue_new == "Blood/Immune", paste0(tissue,".",cell_ontology),cell_ontology)) 

#plot the four heatmaps
col_fun <- colorRamp2(seq(-0.3,1,0.1), viridis::viridis(14))
ht_opt(TITLE_PADDING=unit(5,"mm"))

facs_mean_hm <- Heatmap(facs_mean_cor, name = "correlation",cluster_rows = T,cluster_columns = T,show_column_dend = F, show_row_dend = F,rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
                        column_split = facs_ct_meta$tissue_new[match(colnames(facs_mean_cor), facs_ct_meta$cell_type)],
                        column_labels =facs_ct_meta$ontology_new[match(colnames(facs_mean_cor),facs_ct_meta$cell_type)],
                        row_split = facs_ct_meta$tissue_new[match(colnames(facs_mean_cor), facs_ct_meta$cell_type)],
                        row_labels =facs_ct_meta$ontology_new[match(colnames(facs_mean_cor),facs_ct_meta$cell_type)],
                        row_names_side = "right",
                        column_names_rot = 55,
                        column_title_rot = 30,
                        row_title_rot = 60,
                        column_title_gp = gpar(fontsize=9),
                        row_title_gp = gpar(fontsize = 9),
                        row_names_gp = gpar(fontsize=8),
                        column_names_gp = gpar(fontsize=8))

droplet_mean_hm <- Heatmap(droplet_mean_cor, name = "correlation",cluster_rows = T,cluster_columns = T,show_column_dend = F, show_row_dend = F,rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
                        column_split = droplet_ct_meta$tissue_new[match(colnames(droplet_mean_cor), droplet_ct_meta$cell_type)],
                        column_labels = droplet_ct_meta$ontology_new[match(colnames(droplet_mean_cor),droplet_ct_meta$cell_type)],
                        row_split = droplet_ct_meta$tissue_new[match(colnames(droplet_mean_cor), droplet_ct_meta$cell_type)],
                        row_labels = droplet_ct_meta$ontology_new[match(colnames(droplet_mean_cor),droplet_ct_meta$cell_type)],
                        row_names_side = "right",
                        column_names_rot = 55,
                        column_title_rot = 30,
                        row_title_rot = 60,
                        column_title_gp = gpar(fontsize=9),
                        row_title_gp = gpar(fontsize = 9),
                        row_names_gp = gpar(fontsize=8),
                        column_names_gp = gpar(fontsize=8))

facs_tissue_factor = factor(facs_ct_meta$tissue_new[match(colnames(facs_mean_cor), facs_ct_meta$cell_type)], levels = names(row_order(facs_mean_hm)))
facs_sscore_hm <- Heatmap(facs_sscore_cor, name = "correlation",cluster_rows = F,cluster_columns = F,show_column_dend = F, show_row_dend = F,rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
                        column_split = facs_tissue_factor,
                        column_labels =facs_ct_meta$ontology_new[match(colnames(facs_mean_cor),facs_ct_meta$cell_type)],
                        row_split =  facs_tissue_factor,
                        row_labels =facs_ct_meta$ontology_new[match(colnames(facs_mean_cor),facs_ct_meta$cell_type)],
                        row_names_side = "right",
                        column_names_rot = 55,
                        column_title_rot = 30,
                        row_title_rot = 60,
                        row_order = unlist(row_order(facs_mean_hm)), 
                        column_order = unlist(column_order(facs_mean_hm)),
                        column_title_gp = gpar(fontsize=9),
                        row_title_gp = gpar(fontsize = 9),
                        row_names_gp = gpar(fontsize=8),
                        column_names_gp = gpar(fontsize=8))

droplet_tissue_factor = factor(droplet_ct_meta$tissue_new[match(colnames(droplet_mean_cor), droplet_ct_meta$cell_type)], levels = names(row_order(droplet_mean_hm)))
droplet_sscore_hm <- Heatmap(droplet_sscore_cor, name = "correlation",cluster_rows = F,cluster_columns = F,show_column_dend = F, show_row_dend = F,rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
                           column_split = droplet_tissue_factor,
                           column_labels = droplet_ct_meta$ontology_new[match(colnames(droplet_mean_cor),droplet_ct_meta$cell_type)],
                           row_split = droplet_tissue_factor,
                           row_labels = droplet_ct_meta$ontology_new[match(colnames(droplet_mean_cor),droplet_ct_meta$cell_type)],
                           row_names_side = "right",
                           column_names_rot = 55,
                           column_title_rot = 30,
                           row_title_rot = 60,
                           row_order = unlist(row_order(droplet_mean_hm)), 
                           column_order = unlist(column_order(droplet_mean_hm)),
                           column_title_gp = gpar(fontsize=9),
                           row_title_gp = gpar(fontsize = 9),
                           row_names_gp = gpar(fontsize=8),
                           column_names_gp = gpar(fontsize=8))

#### maximize cell type specificity and expressed cell types ####
max_specificity <- rowMax(facs_sscore) %>% set_names(rownames(facs_sscore))
num_ct_expression <- rowSums(as.matrix(t(facs_mean))>0 ) %>% .[match(names(max_specificity), names(.))]

#calculate entropy
facs_discrete_mat <- as.matrix(assay(facs_obj_sce, "logcounts"))[match(names(max_specificity), rownames(facs_obj_sce)), ] 
facs_gene_entropy <- map(rownames(facs_discrete_mat), ~infotheo::discretize(facs_discrete_mat[.x, ])) %>%
  map(~infotheo::entropy(.x$X))

facs_gene_ct_entropy <- as.matrix(t(facs_mean))[match(names(max_specificity), colnames(facs_mean)), ] 
facs_gene_ct_entropy <-  map(rownames(facs_gene_ct_entropy), ~infotheo::discretize(facs_gene_ct_entropy[.x, ])) %>%
  map(~infotheo::entropy(.x$X))

##### Simulated data results####
#load simulated score
parameter_df <- read.table(here("data", "expr", "score_sim", "parameter_df.txt"), header=T) %>%
  as_tibble()

all_res <- map(parameter_df$output_header, ~ {if (file.exists(paste0(.x, ".score_df.txt" ))) read.csv(paste0(.x, ".score_df.txt" )) else NA}) %>%
  map(~dplyr::rename(.x, "gene" = "X")) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = .y)) %>%
  map2(parameter_df$rep, ~mutate(.x, rep = .y)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  group_by(gene, effect_size) %>% 
  summarise(mean_seismic_score = mean(seismic_score), sd_seismic_score = sd(seismic_score), 
            mean_mean_value = mean(mean_value), sd_mean_value = sd(mean_value),
            mean_spc_score = mean(spc_score), sd_spc_score = sd(spc_score),
            mean_de_score = mean(de_score), sd_de_score = sd(de_score) ) %>% 
  ungroup() %>%
  mutate(effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = effect_size)))

#load gene mapping
load(here("data","expr","null_sim","expr_rda_rs","expr_ds_1.rda"))
load(here("data","expr","Tabula_muris","facs_clean.rda") )
facs_obj <- facs_obj_sce[,which(!is.na(facs_obj_sce$cell_ontology_class))]

set.seed(99)
cell_idx <- sample(1:ncol(facs_obj), size=10000) #select 10000 cells
gene_idx <- sample(1:nrow(facs_obj)) #gene index
gene_mapping = tibble(original_name = rownames(facs_obj), new_names = rownames(facs_obj)[gene_idx])

#map to actual gene names
all_res <- all_res %>%
  left_join(gene_mapping, by = c("gene" = "new_names")) %>%
  select(-gene) %>%
  dplyr::rename("original_name" = "gene")

#plot heatmap
all_res_df <- all_res %>% 
  select(-contains("sd")) %>%
  pivot_longer( cols = starts_with("mean_"),names_to = "metric",values_to = "value", names_prefix = "mean_") %>%
  arrange(effect_size) 

quantile_data <- all_res_df %>%
  filter(effect_size == 0) %>%
  group_by(metric) %>%
  mutate(quantile = ntile(value, 5) ) %>%
  mutate(quantile = paste0("Q", quantile)) %>% 
  ungroup() %>%
  select(gene, metric, quantile)

all_res_df <- all_res_df %>%
  left_join(quantile_data, by = c("gene", "metric") ) %>%
  mutate(metric = recode( metric, "spc_score" = "specificity", "seismic_score" = "seismic score", "mean_value" = "mean expression", "de_score" = "deScore")) %>%
  mutate(metric = factor(metric, levels = c("seismic score", "deScore", "mean expression", "specificity")))

ggplot(all_res_df, aes(x = effect_size, y = value, group = gene, color = metric)) +
  geom_line(alpha = 0.5) +
  ggsci::scale_color_npg()+
  facet_wrap(~  metric + quantile, scales = "free", ncol = 5)  +
  theme_classic() +
  xlab("normalized effect size")

#selected plots
selected_groups <- all_res %>% 
  filter(gene %in%  sample(unique(gene), 10)) %>%
  mutate(effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = effect_size))) %>%
  pivot_longer( cols = -c(effect_size, gene),
                names_to = c("stat", "score_type"),
                names_pattern = "^(mean|sd)_(.*)$") %>%
  pivot_wider(names_from  = stat, values_from = value) 

selected_groups <- selected_groups %>%
  mutate(mean_max = mean + sd, mean_min = mean - sd)

ggplot(selected_groups, aes(x = effect_size, y = mean, group = gene, ymin = mean_min, ymax=mean_max,)) +
  geom_point() +
  geom_line() +
  geom_ribbon(alpha=0.5) + 
  facet_wrap(~ gene + score_type, scales = "free", ncol = 4) 

#load simulated scores
parameter_df <- read.table(here("data", "expr", "score_sim", "parameter_df.txt"), header=T) %>%
  as_tibble()

seismic_score_res <- paste0(parameter_df$output_header, ".seismic_score.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~ {if (!is.data.frame(.x)) NA else pivot_longer(.x, cols = !gene, names_to = "cell_type")}) %>% 
  map2(parameter_df$es_name, ~ {if (!is_tibble(.x)) NA else mutate(.x, effect_size = .y)}) %>%
  keep(~is_tibble(.)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>% 
  group_by(effect_size, gene, cell_type) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  distinct()

seismic_p_res <- paste0(parameter_df$output_header, ".seismic_p.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~ {if (!is.data.frame(.x)) NA else pivot_longer(.x, cols = !gene, names_to = "cell_type")}) %>% 
  map2(parameter_df$es_name, ~ {if (!is_tibble(.x)) NA else mutate(.x, effect_size = .y)}) %>%
  keep(~is_tibble(.)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>% 
  group_by(effect_size, gene, cell_type) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  distinct()

seismic_r_res <- paste0(parameter_df$output_header, ".seismic_r.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~ {if (!is.data.frame(.x)) NA else pivot_longer(.x, cols = !gene, names_to = "cell_type")}) %>% 
  map2(parameter_df$es_name, ~ {if (!is_tibble(.x)) NA else mutate(.x, effect_size = .y)}) %>%
  keep(~is_tibble(.)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>% 
  group_by(effect_size, gene, cell_type) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  distinct()

descore_res <- paste0(parameter_df$output_header, ".de_score.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~ {if (!is.data.frame(.x)) NA else pivot_longer(.x, cols = !gene, names_to = "cell_type")}) %>% 
  map2(parameter_df$es_name, ~ {if (!is_tibble(.x)) NA else mutate(.x, effect_size = .y)}) %>%
  keep(~is_tibble(.)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>% 
  group_by(effect_size, gene, cell_type) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  distinct()

mean_res <- paste0(parameter_df$output_header, ".mean_val.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~ {if (!is.data.frame(.x)) NA else pivot_longer(.x, cols = !gene, names_to = "cell_type")}) %>% 
  map2(parameter_df$es_name, ~ {if (!is_tibble(.x)) NA else mutate(.x, effect_size = .y)}) %>%
  keep(~is_tibble(.)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>% 
  group_by(effect_size, gene, cell_type) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  distinct()

spc_score_res <- paste0(parameter_df$output_header, ".spc_score.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~ {if (!is.data.frame(.x)) NA else pivot_longer(.x, cols = !gene, names_to = "cell_type")}) %>% 
  map2(parameter_df$es_name, ~ {if (!is_tibble(.x)) NA else mutate(.x, effect_size = .y)}) %>%
  keep(~is_tibble(.)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>% 
  group_by(effect_size, gene, cell_type) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  distinct()

seismic_score_res <- seismic_score_res %>% mutate(score_type = "seismic specificity score")
seismic_p_res <- seismic_p_res %>%  mutate(score_type = "seismic p score")
seismic_r_res <- seismic_r_res %>%  mutate(score_type = "seismic r score")
descore_res <- descore_res %>% filter(value < 500) %>% mutate(score_type = "descore")
mean_res <- mean_res %>% mutate(score_type = "mean expression")
spc_score_res <- spc_score_res %>% mutate(score_type = "Bryois score")

all_res_score <- rbind(seismic_score_res, seismic_p_res, seismic_r_res, descore_res, mean_res, spc_score_res) %>%
  mutate(is_target = ifelse(cell_type == "target_cell_type", "target cell type", "other cell types")) %>%
  mutate(effect_size = as.numeric(gsub(pattern = "es_", x = effect_size, replacement = "")))

ggplot(all_res_score, aes(x = factor(effect_size), y = value, fill = is_target)) +
  geom_boxplot(position = position_dodge2(width = 0.2, preserve = "single"),outlier.shape = NA) +
  facet_wrap(~ score_type, scales = "free_y") +
  theme_classic() +
  ggh4x::facetted_pos_scales(
    y = list(
      score_type == "Bryois score" ~ scale_y_continuous(limits = c(0,0.5)),
      score_type == "descore" ~ scale_y_continuous(limits = c(0,10)),
      score_type == "mean expression" ~ scale_y_continuous(limits = c(0,40)),
      score_type == "seismic p score" ~ scale_y_continuous(limits = c(0,1)),
      score_type == "seismic r score" ~ scale_y_continuous(limits = c(0,1)),
      score_type == "seismic specificity score" ~ scale_y_continuous(limits = c(0,0.1))
    )
  )

##new
parameter_df <- read.table(here("data", "expr", "score_causal", "parameter_df.txt"), header=T) %>%
  as_tibble()

#load replace_mat 
target_gene_list <- parameter_df$perturbed_mat_file %>%
  map(~read.table(.x)) %>%
  map(~rownames(.x))
  
seismic_score_res <- paste0(parameter_df$output_header, ".seismic_score.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  #map(~select(.x, gene, target_cell_type)) %>%
  map(~pivot_longer(.x, cols = !gene, names_to = "cell_type")) %>%
  map2(target_gene_list, ~mutate(.x, gene_type = ifelse(gene %in% .y, "target genes", "others"))) %>%
  map(~drop_na(.x)) 

de_score_res <- paste0(parameter_df$output_header, ".de_score.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~pivot_longer(.x, cols = !gene, names_to = "cell_type")) %>%
  map2(target_gene_list, ~mutate(.x, gene_type = ifelse(gene %in% .y, "target genes", "others"))) %>%
  map(~drop_na(.x)) 

mean_res <- paste0(parameter_df$output_header, ".mean_val.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~pivot_longer(.x, cols = !gene, names_to = "cell_type")) %>%
  map2(target_gene_list, ~mutate(.x, gene_type = ifelse(gene %in% .y, "target genes", "others"))) %>%
  map(~drop_na(.x)) 


#get the test
seismic_score_mi <- seismic_score_res %>%
  #map(~infotheo::mutinformation( infotheo::discretize(.x$target_cell_type), .x$gene_type)) 
  map(~ group_by(.x, cell_type)) %>%
  map(~ mutate(.x, discrete_value = infotheo::discretize(value))) %>%
  map(~ summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type)))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  map2(parameter_df$gr_name, ~mutate(.x, gr_name = .y)) %>% 
  map2(parameter_df$gs_name, ~mutate(.x, gs_name = .y)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "seismic specificity")

de_score_mi <- de_score_res %>% 
  # map(~infotheo::mutinformation( infotheo::discretize(.x$target_cell_type), .x$gene_type)) 
  map(~ group_by(.x, cell_type)) %>%
  map(~ mutate(.x, discrete_value = infotheo::discretize(value))) %>%
  map(~ summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type)))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  map2(parameter_df$gr_name, ~mutate(.x, gr_name = .y)) %>% 
  map2(parameter_df$gs_name, ~mutate(.x, gs_name = .y)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "DE score")

mean_mi <- mean_res %>% 
 # map(~infotheo::mutinformation( infotheo::discretize(.x$target_cell_type), .x$gene_type)) 
  map(~ group_by(.x, cell_type)) %>%
  map(~ mutate(.x, discrete_value = infotheo::discretize(value))) %>%
  map(~ summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type)))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  map2(parameter_df$gr_name, ~mutate(.x, gr_name = .y)) %>% 
  map2(parameter_df$gs_name, ~mutate(.x, gs_name = .y)) %>% 
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "mean expression")

mi_res <- rbind(seismic_score_mi, de_score_mi, mean_mi) %>%
  group_by(effect_size, gr_name, gs_name, score_type) %>%
  mutate(nmi = nmi - (sum(nmi)-nmi)/(length(nmi) - 1)) %>%
  ungroup() %>%
  filter(cell_type == "target_cell_type")
  

ggplot(mi_res,  aes(x = factor(effect_size), y = nmi, fill = score_type)) +
  geom_boxplot(position = position_dodge2(width = 0.1, preserve = "single"), alpha = 0.8) +
  #facet_wrap(~score_type,ncol = 1) +
  theme_classic() +
  ggsci::scale_fill_npg()


#new mi
seismic_score_mi <- seismic_score_res %>%
  #map(~infotheo::mutinformation( infotheo::discretize(.x$target_cell_type), .x$gene_type)) 
  map(~ group_by(.x, cell_type)) %>%
  map(~ mutate(.x, discrete_value = infotheo::discretize(value))) %>%
  map(~ summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type)))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "seismic specificity")

de_score_mi <- de_score_res %>% 
  # map(~infotheo::mutinformation( infotheo::discretize(.x$target_cell_type), .x$gene_type)) 
  map(~ group_by(.x, cell_type)) %>%
  map(~ mutate(.x, discrete_value = infotheo::discretize(value))) %>%
  map(~ summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type)))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "DE score")

mean_mi <- mean_res %>% 
  # map(~infotheo::mutinformation( infotheo::discretize(.x$target_cell_type), .x$gene_type)) 
  map(~ group_by(.x, cell_type)) %>%
  map(~ mutate(.x, discrete_value = infotheo::discretize(value))) %>%
  map(~ summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type)))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "mean expression")

mi_res <- rbind(seismic_score_mi, de_score_mi, mean_mi) %>%
  mutate(perturbed_cell_type = ifelse(cell_type == "target_cell_type", T, F))

ggplot(mi_res,  aes(x = factor(effect_size), y = nmi, fill = perturbed_cell_type)) +
  geom_boxplot(position = position_dodge2(width = 0.1, preserve = "single"), alpha = 0.8) +
  facet_wrap(~score_type,ncol = 1) +
  theme_classic() +
  ggsci::scale_color_npg()


#do the test
library(fgsea)
seismic_score_p <- seismic_score_res %>%
  map(~setNames(.x$target_cell_type, .x$gene)) %>%
  map2(seismic_score_res, ~fgsea(list(target_genes = .y$gene[.y$gene_type != "others"]), .x))

seismic_r_p <- seismic_r_res %>%
  map(~setNames(.x$target_cell_type, .x$gene)) %>%
  map2(seismic_r_res, ~fgsea(list(target_genes = .y$gene[.y$gene_type != "others"]), .x))

seismic_p_p <- seismic_p_res %>%
  map(~setNames(.x$target_cell_type, .x$gene)) %>%
  map2(seismic_p_res, ~fgsea(list(target_genes = .y$gene[.y$gene_type != "others"]), .x))

spc_score_p <- spc_score_res %>%
  map(~setNames(.x$target_cell_type, .x$gene)) %>%
  map2(spc_score_res, ~fgsea(list(target_genes = .y$gene[.y$gene_type != "others"]), .x))

mean_p <- mean_res %>% 
  map(~setNames(.x$target_cell_type, .x$gene)) %>%
  map2(mean_res, ~fgsea(list(target_genes = .y$gene[.y$gene_type != "others"]), .x))

p_res <- parameter_df %>% 
  select(es_name) %>%
  mutate(seismic_score = unlist(seismic_score_p %>% map(~.x$pval)) , 
         seismic_r = unlist(seismic_r_p %>% map(~.x$pval)),  
         seismic_p = unlist(seismic_p_p %>% map(~.x$pval)), 
         spc_score = unlist(spc_score_p %>% map(~.x$pval)), 
         mean_exp = unlist(mean_p%>% map(~.x$pval))) %>%
  group_by(es_name) %>%
  mutate_if(is.numeric, ~p.adjust(.x, method = "fdr")) %>%
  mutate(effect_size = as.numeric(gsub(pattern = "es_", x = es_name, replacement = ""))) %>%
  ungroup() %>%
  select(-es_name) %>%
  pivot_longer(cols = !effect_size, names_to = "score_type", values_to = "value")%>% 
  group_by(effect_size, score_type) %>% summarise(ratio = length(which(value<=0.05))/ length(value))

ggplot(p_res %>% ungroup(),  aes(x = factor(effect_size), y = ratio, group = score_type, color = score_type)) +
  geom_line() +
  theme_classic()

###all paramter df 
parameter_df <- c("causal_overlap" = here("data", "expr", "score_causal", "multi_ct", "causal_overlap", "parameter_df.txt"),
                  "ct_overlap" = here("data", "expr", "score_causal", "multi_ct", "ct_overlap", "parameter_df.txt"),
                  "no_overlap" = here("data", "expr", "score_causal", "multi_ct", "no_overlap", "parameter_df.txt"),
                  "random" = here("data", "expr", "score_causal", "multi_ct", "random", "parameter_df.txt")) %>%
  map(~read.table(.x, header=T)) %>% 
  map2(names(.), ~mutate(.x, simulation_type = .y)) %>%
  purrr::reduce(~rbind(.x, .y))

#load replace_mat  (target gene list)
# target_gene_list <- parameter_df$perturbed_mat_file %>%
#   map(~read.table(.x)) %>%
#   map(~rownames(.x))

all_gene_anno <-  parameter_df$gene_anno_file %>%
  map(~read.table(.x, header = T)) %>% 
  map(~ {data_tbl <- .x; map(c(1,2,3), ~select(data_tbl, any_of(c("mmu_symbol", "is_ct_specific", paste0("is_ct_specific_", .x), "is_causal", paste0("is_causal_", .x)))))}) %>%
  map(~map(.x, ~filter(.x, .[[2]] | .[[3]]))) %>%
  map(~map(.x, ~select(.x, mmu_symbol))) %>%
  map(~map2(.x, c(1,2,3), ~mutate(.x, cell_type = paste0("target_cell_type_", .y)))) %>%
  map(~purrr::reduce(.x, ~rbind(.x, .y))) %>%
  map(~rename(.x, "gene" = "mmu_symbol"))
  
#load results
seismic_score_res <- paste0(parameter_df$output_header, ".seismic_score.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~pivot_longer(.x, cols = !gene, names_to = "cell_type")) %>%
  map2(all_gene_anno, ~mutate(.x, gene_type = ifelse(paste0(gene, cell_type) %in% paste0(.y[["gene"]], .y[["cell_type"]]), "target genes", "others"))) %>%
  map(~drop_na(.x)) 

de_score_res <- paste0(parameter_df$output_header, ".de_score.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~ {if (is.data.frame(.x)) pivot_longer(.x, cols = !gene, names_to = "cell_type") else NA}) %>%
  map2(all_gene_anno, ~{if (is_tibble(.x)) mutate(.x, gene_type = ifelse(paste0(gene, cell_type) %in% paste0(.y[["gene"]], .y[["cell_type"]]), "target genes", "others")) else NA}) %>%
  map(~ {if (is_tibble(.x)) drop_na(.x) else NA}) 

mean_res <- paste0(parameter_df$output_header, ".mean_val.txt") %>%
  map(~ {if (file.exists(.x)) read.table(.x, header = T, sep = "\t") else NA}) %>%
  map(~ pivot_longer(.x, cols = !gene, names_to = "cell_type")) %>%
  map2(all_gene_anno, ~mutate(.x, gene_type = ifelse(paste0(gene, cell_type) %in% paste0(.y[["gene"]], .y[["cell_type"]]), "target genes", "others"))) %>%
  map(~drop_na(.x)) 

#calculate mi
seismic_score_mi <- seismic_score_res %>%
  map(~filter(.x, grepl(pattern = "target_cell_type", x = cell_type))) %>%
  map(~group_by(.x, cell_type)) %>%
  map(~mutate(.x, discrete_value = infotheo::discretize(value))) %>%
  map(~summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type)))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  map2(parameter_df$gs_name, ~mutate(.x, gs_name = .y)) %>% 
  map2(parameter_df$simulation_type, ~mutate(.x, simulation_type = .y)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "seismic specificity")

de_score_mi <- de_score_res %>% 
  map(~ {if (is_tibble(.x)) filter(.x, grepl(pattern = "target_cell_type", x = cell_type)) NA}) %>%
  map(~ {if (is_tibble(.x)) group_by(.x, cell_type) else NA})%>%
  map(~ {if (is_tibble(.x)) mutate(.x, discrete_value = infotheo::discretize(value)) else NA}) %>%
  map(~ {if (is_tibble(.x)) summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type))) else NA}) %>%
  map2(parameter_df$es_name, ~{if (is_tibble(.x)) mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y))) else NA}) %>% 
  map2(parameter_df$gs_name, ~{if (is_tibble(.x)) mutate(.x, gs_name = .y) else NA}) %>% 
  map2(parameter_df$simulation_type, ~{if (is_tibble(.x)) mutate(.x, simulation_type = .y) else NA}) %>%
  keep(!is.na(.)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "DE score")

mean_mi <- mean_res %>% 
  map(~filter(.x, grepl(pattern = "target_cell_type", x = cell_type))) %>%
  map(~group_by(.x, cell_type)) %>%
  map(~mutate(.x, discrete_value = infotheo::discretize(value))) %>%
  map(~summarise(.x, nmi = infotheo::mutinformation(discrete_value, gene_type)/ min(infotheo::entropy(discrete_value), infotheo::entropy(gene_type)))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  map2(parameter_df$gs_name, ~mutate(.x, gs_name = .y)) %>% 
  map2(parameter_df$simulation_type, ~mutate(.x, simulation_type = .y)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "mean expression")

mi_res <- rbind(seismic_score_mi, de_score_mi, mean_mi) %>%
  mutate(perturbed_cell_type = ifelse(grepl(pattern = "target_cell_type", x = cell_type), T, F))

ggplot(mi_res,  aes(x = factor(effect_size), y = nmi, fill = perturbed_cell_type)) +
  geom_boxplot(position = position_dodge2(width = 0.1, preserve = "single"), alpha = 0.8) +
  facet_wrap(~score_type + simulation_type,ncol = 4) +
  theme_classic() +
  ggsci::scale_color_npg()

#plot null
mean_res_null <- mean_res %>%
  map(~filter(.x, gene_type != "target genes", grepl(pattern = "target_cell_type", x = cell_type))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  map2(parameter_df$gs_name, ~mutate(.x, gs_name = .y)) %>% 
  map2(parameter_df$simulation_type, ~mutate(.x, simulation_type = .y)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "mean expression")

de_score_res_null <- de_score_res %>% 
  map(~ {if (is_tibble(.x)) filter(.x, gene_type != "target genes", grepl(pattern = "target_cell_type", x = cell_type)) else NA}) %>%
  map2(parameter_df$es_name, ~ {if (is_tibble(.x)) mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y))) else NA}) %>% 
  map2(parameter_df$gs_name, ~ {if (is_tibble(.x)) mutate(.x, gs_name = .y) else NA}) %>% 
  map2(parameter_df$simulation_type, ~ {if (is_tibble(.x)) mutate(.x, simulation_type = .y) else NA}) %>%
  keep(!is.na(.)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "DE score")

seismic_score_null <- seismic_score_res %>% 
  map(~filter(.x, gene_type != "target genes", grepl(pattern = "target_cell_type", x = cell_type))) %>%
  map2(parameter_df$es_name, ~mutate(.x, effect_size = as.numeric(gsub(pattern = "es_", replacement = "", x = .y)))) %>% 
  map2(parameter_df$gs_name, ~mutate(.x, gs_name = .y)) %>% 
  map2(parameter_df$simulation_type, ~mutate(.x, simulation_type = .y)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(score_type = "seismic specificity")

null_res <- rbind(mean_res_null, de_score_res_null, seismic_score_null) 

ggplot(null_res %>% slice_sample(n=100000),  aes(x = factor(effect_size), y = value)) +
  geom_boxplot(position = position_dodge2(width = 0.1, preserve = "single"), alpha = 0.8) +
  facet_wrap(~ score_type + simulation_type , nrow= 3 , scales = "free") +
  theme_classic() 


####load simulation score (with different perturbation ratio)
parameter_df <- here("data", "expr", "score_robustness", "target_sample_5.0", "new_parameter_df") %>%
  read.table(header=T) %>%
  as_tibble()

target_gene_list <- parameter_df$perturbed_mat_file %>%
  map(~read.table(.x)) %>%
  map(~rownames(.x))

all_results <- parameter_df$output_header %>%
  map(~ {file_header = .x; map(1:10, ~read.table(paste0(file_header, ".sample_", .x, ".score_df.txt"), header = T))})

#load calculate specificity score of previous 


