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
if (!require("seismicGWAS")){
  if (!requireNamespace("devtools", quietly = TRUE)){
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

##### plot marker distribution in pancreatic cell types#####
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
  map(~{if (length(.x) == 5) .x else c(.x, list(empty_plot)[rep(1, 5 - length(.x))])}) %>%
  purrr::reduce(~c(.x, .y))

final_grid <- final_grid <- cowplot::plot_grid(plotlist = pancreas_sscore_rank_plot_list,ncol = 4, nrow = 5,align = "hv",byrow = F)

##### plot violin plot to compare between housekeeping genes and other target genes #####
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
