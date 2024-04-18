#generate plots for results of the Tabula Muris FACS dataset
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

#encode trait mapping
neuropsy_disesaes = c("insominia","BD","MDD","Nrt_new","Scz_new")
immune_diseases = c("AID","Hypothyroidism","RA","ibd","cd","uc")
others = c("Smoking","BMI","College_edu","RBC","Lymphocyte_count","Monocyte_count","T1D","T2D_2","Cardiovas","SBP","AF","CKD","glucose_q","HDL_q","LDL_q","TG_q")
trait_meta = tibble(trait_names = c(neuropsy_disesaes, immune_diseases, others),
                    official_names = c("Insomnia","Bipolar disorder","Depression","Neuroticism","Schizophrenia",
                                        "Autoimmune diseases","Hypothroidism","Rheumatoid arthritis","Inflammatory bowel disease","Crohn's disease","Ulcerative colitis",
                                       "Smoking","BMI","College education","Erythrocyte count","Lymphocyte count","Monocyte count","Type I diabetes","Type II diabetes","Cardiovascular diseases",
                                       "Systolic blood pressure","Atrial fibrillation","Chronic kidney disease","Glucose level","HDL level","LDL level","Triglyceride level"),
                    trait_type = factor(c(rep("neuropsy",length(neuropsy_disesaes)), rep("immune",length(immune_diseases)), rep("others",length(others))), levels=c("neuropsy","immune","others")) )

#load results - facs
facs_res = read.table( here("results","Tabula_muris","FACS","seismic","facs_res.txt"), header = T, sep = "\t") %>% 
  as_tibble() %>%
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))  %>%
  mutate(cell_type = factor(cell_type, levels = sort(cell_type))) %>%
  arrange(cell_type)

facs_scdrs_res = list.files(here("results","Tabula_muris","FACS","scDRS"),pattern = "*.scdrs_group.cluster_name",full.names = T) %>% 
  set_names(str_extract(string = ., pattern = "(?<=/)[^/]+$") %>% gsub(pattern = ".scdrs_group.cluster_name", replacement = "", x=.)) %>%
  map(~read.table(.x, header=T, sep="\t") %>% as_tibble()) %>%
  map(~mutate(.x, P_val = pnorm(assoc_mcz, lower.tail = F))) %>%
  map(~select(.x, any_of(c("group","P_val")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>% 
  purrr::reduce(~left_join(.x, .y, by = "cell_type")) %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))%>%
  mutate(cell_type = factor(cell_type, levels = levels(facs_res$cell_type))) %>%
  arrange(cell_type)

facs_fuma_res = list.files(here("results","Tabula_muris","FACS","FUMA"), pattern = "*.gsa.out", full.names = T) %>%
  set_names(str_extract(string = . , pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.[A-Za-z.]+$", replacement = "", x= .)) %>%
  map(~read.table(.x, header = T) %>% as_tibble()) %>%
  map(~select(.x, any_of(c("VARIABLE", "P")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by = "cell_type"))  %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) 

facs_magma_res = list.files(here("results","Tabula_muris","FACS","S-MAGMA"), pattern = "*.gsa.out", full.names = T) %>%
  set_names(str_extract(string = . , pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.[A-Za-z.]+$", replacement = "", x= .)) %>%
  map(~read.table(.x, header = T) %>% as_tibble()) %>%
  map(~select(.x, any_of(c("VARIABLE", "P")))) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by = "cell_type"))  %>% 
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) 

#load auxillary files of fuma/s-magma and harmonise cell type names
facs_ct_meta = tibble(cell_type = as.character(facs_res$cell_type)) %>%
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
  mutate(tissue_new = ifelse(grepl(pattern = "B cell|T cell|immune cells|erythrocyte|macrophage|microglia|myeloid|leukocyte|natural killer|antigen|NK|monocyte|lymphocyte|blood|thymocyte", x=cell_ontology), "Blood/Immune",tissue_new)) %>%
  mutate(ontology_new = ifelse(tissue_new == "Blood/Immune", paste0(tissue,".",cell_ontology),cell_ontology)) 

#harmonise cell types for FUMA/S-MAGMA output
facs_fuma_ct_mapping = read.table(here("data","expr","Tabula_muris","tm_facs.fuma.aux.txt"), header = T, sep = "\t") %>% as_tibble()
facs_magma_ct_mapping = read.table(here("data","expr","Tabula_muris","tm_facs.magma.aux.txt"), header = T, sep = "\t") %>% as_tibble()

facs_fuma_res = facs_fuma_ct_mapping %>% 
  left_join(facs_fuma_res, by = c("encoded_name" = "cell_type")) %>%
  select( -encoded_name )  %>% 
  mutate(cell_type = factor(cell_type, levels = levels(facs_res$cell_type))) %>%
  arrange(cell_type)

facs_magma_res = facs_magma_ct_mapping %>% 
  left_join(facs_magma_res, by = c("encoded_name" = "cell_type")) %>%
  select( -encoded_name ) %>% 
  mutate(cell_type = factor(cell_type, levels = levels(facs_res$cell_type))) %>%
  arrange(cell_type)


##### plot venn diagram ######
#figure 3 B
library("ggvenn")
facs_seismic_pair = facs_res %>%
  mutate(across(where(is.double) , ~p.adjust(.x, method ="fdr"))) %>%
  pivot_longer(!cell_type, names_to = "trait", values_to = "fdr") %>%
  filter(fdr<=0.05) %>% mutate(pair = paste0(trait,".",cell_type)) %>% 
  pull(pair)

facs_scdrs_pair = facs_scdrs_res %>%
  mutate(across(where(is.double) , ~p.adjust(.x, method ="fdr"))) %>%
  pivot_longer(!cell_type, names_to = "trait", values_to = "fdr") %>% 
  filter(fdr<=0.05) %>% mutate(pair = paste0(trait,".",cell_type)) %>% 
  pull(pair)  

facs_fuma_pair = facs_fuma_res %>%
  mutate(across(where(is.double) , ~p.adjust(.x, method ="fdr"))) %>%
  pivot_longer(!cell_type, names_to = "trait", values_to = "fdr") %>% 
  filter(fdr<=0.05) %>% 
  mutate(pair = paste0(trait,".",cell_type)) %>% 
  pull(pair)

facs_magma_pair = facs_magma_res %>%
  mutate(across(where(is.double) , ~p.adjust(.x, method ="fdr"))) %>%
  pivot_longer(!cell_type, names_to = "trait", values_to = "fdr") %>% 
  filter(fdr<=0.05) %>% 
  mutate(pair = paste0(trait,".",cell_type)) %>% 
  pull(pair)

ordered_colors = c(color_mapping_vec["S-MAGMA"], color_mapping_vec["FUMA"], color_mapping_vec["scDRS"], color_mapping_vec["seismic"])
facs_venn = list(`S-MAGMA` = facs_magma_pair, FUMA = facs_fuma_pair, scDRS = facs_scdrs_pair, seismic = facs_seismic_pair) 
ggvenn(facs_venn,text_size = 6,set_name_size=6,show_percentage = F,fill_color = ordered_colors %>% set_names(NULL)  ,fill_alpha = 0.4) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-2, 1.5))


#### Plot heatmap #####
#figure 3A  + sup fig 3-5
#define paramters
library(ComplexHeatmap)
library(circlize)
ht_opt(TITLE_PADDING=unit(5,"mm"))
col_fun = colorRamp2(c(-0.5, -log10(0.05), 4,10), viridis::viridis(4))

#seismic heatmap
facs_fdr_mat = facs_res %>%
  mutate(across(where(is.double), ~p.adjust(.x, method ="fdr"))) %>%
  select(all_of(c("cell_type",neuropsy_disesaes,immune_diseases,others))) %>%
  select(where(is.double)) %>%
  as.matrix() %>%
  t() %>%
  set_colnames(facs_ct_meta$cell_type[match(facs_res$cell_type,facs_ct_meta$cell_type)]) %>%
  set_rownames(trait_meta$official_names[match(rownames(.), trait_meta$trait_names)])

Heatmap(-log10(facs_fdr_mat), name = "-log10(FDR)",cluster_rows = F,cluster_columns = F,show_column_dend = F, show_row_dend = F,rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
        column_split =facs_ct_meta$tissue_new[match(colnames(facs_fdr_mat),facs_ct_meta$cell_type)],
        column_labels =facs_ct_meta$ontology_new[match(colnames(facs_fdr_mat),facs_ct_meta$cell_type)],
        row_split = trait_meta$trait_type[match(rownames(facs_fdr_mat), trait_meta$official_names)] %>% factor(levels = c("neuropsy","immune","others")),
        row_names_side = "left",
        column_title_gp = gpar(fontsize=11),
        row_names_gp = gpar(fontsize=9),
        column_names_gp = gpar(fontsize=8),
        column_names_rot = 55,
        column_title_rot = 30,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(facs_fdr_mat[i, j] < 1e-2) {
            grid.text("**", x, y, gp = gpar(fontsize = 8))
          } else if(facs_fdr_mat[i, j] <= 0.05) {
            grid.text("*", x, y, gp = gpar(fontsize = 8))
          }
        })

# scDRS heatmap
facs_scdrs_mat = facs_scdrs_res %>%
    mutate(across(where(is.double), ~p.adjust(.x, method ="fdr"))) %>%
    select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) %>%
    select(where(is.double)) %>%
    as.matrix() %>%
    t() %>%
    set_colnames(facs_ct_meta$cell_type[match(facs_scdrs_res$cell_type, facs_ct_meta$cell_type)]) %>%
    set_rownames(trait_meta$official_names[match(rownames(.), trait_meta$trait_names)]) %>%
    .[match(rownames(facs_fdr_mat),rownames(.)), match(colnames(facs_fdr_mat), colnames(.))]

Heatmap(-log10(facs_scdrs_mat), name = "-log10(FDR)", cluster_rows = F, cluster_columns = F, show_column_dend = F, show_row_dend = F, rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
        column_split = facs_ct_meta$tissue_new[match(colnames(facs_scdrs_mat), facs_ct_meta$cell_type)],
        column_labels = facs_ct_meta$ontology_new[match(colnames(facs_scdrs_mat), facs_ct_meta$cell_type)],
        row_split = trait_meta$trait_type[match(rownames(facs_scdrs_mat), trait_meta$official_names)] %>% factor(levels = c("neuropsy", "immune", "others")),
        row_names_side = "left",
        column_title_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 55,
        column_title_rot = 30,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (facs_scdrs_mat[i, j] < 1e-2) {
              grid.text("**", x, y, gp = gpar(fontsize = 8))
          } else if (facs_scdrs_mat[i, j] <= 0.05) {
              grid.text("*", x, y, gp = gpar(fontsize = 8))
              }
          })

# FUMA heatmap
facs_fuma_mat = facs_fuma_res %>%
  mutate(across(where(is.double), ~p.adjust(.x, method ="fdr"))) %>%
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) %>%
  select(where(is.double)) %>%
  as.matrix() %>%
  t() %>%
  set_colnames(facs_ct_meta$cell_type[match(facs_fuma_res$cell_type, facs_ct_meta$cell_type)]) %>%
  set_rownames(trait_meta$official_names[match(rownames(.), trait_meta$trait_names)]) %>%
  .[match(rownames(facs_fdr_mat),rownames(.)), match(colnames(facs_fdr_mat), colnames(.))]

Heatmap(-log10(facs_fuma_mat), name = "-log10(FDR)", cluster_rows = F, cluster_columns = F, show_column_dend = F, show_row_dend = F, rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
        column_split = facs_ct_meta$tissue_new[match(colnames(facs_fuma_mat), facs_ct_meta$cell_type)],
        column_labels = facs_ct_meta$ontology_new[match(colnames(facs_fuma_mat), facs_ct_meta$cell_type)],
        row_split = trait_meta$trait_type[match(rownames(facs_fuma_mat), trait_meta$official_names)] %>% factor(levels = c("neuropsy", "immune", "others")),
        row_names_side = "left",
        column_title_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 55,
        column_title_rot = 30,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (facs_fuma_mat[i, j] < 1e-2) {
            grid.text("**", x, y, gp = gpar(fontsize = 8))
          } else if (facs_fuma_mat[i, j] <= 0.05) {
            grid.text("*", x, y, gp = gpar(fontsize = 8))
          }
        })

# S-MAGMA heatmap
facs_magma_mat = facs_magma_res %>%
  mutate(across(where(is.double), ~p.adjust(.x, method ="fdr"))) %>%
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases, others))) %>%
  select(where(is.double)) %>%
  as.matrix() %>%
  t() %>%
  set_colnames(facs_ct_meta$cell_type[match(facs_magma_res$cell_type, facs_ct_meta$cell_type)]) %>%
  set_rownames(trait_meta$official_names[match(rownames(.), trait_meta$trait_names)]) %>%
  .[match(rownames(facs_fdr_mat),rownames(.)), match(colnames(facs_fdr_mat), colnames(.))]

Heatmap(-log10(facs_magma_mat), name = "-log10(FDR)", cluster_rows = F, cluster_columns = F, show_column_dend = F, show_row_dend = F, rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
        column_split = facs_ct_meta$tissue_new[match(colnames(facs_magma_mat), facs_ct_meta$cell_type)],
        column_labels = facs_ct_meta$ontology_new[match(colnames(facs_magma_mat), facs_ct_meta$cell_type)],
        row_split = trait_meta$trait_type[match(rownames(facs_magma_mat), trait_meta$official_names)] %>% factor(levels = c("neuropsy", "immune", "others")),
        row_names_side = "left",
        column_title_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 55,
        column_title_rot = 30,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (facs_magma_mat[i, j] < 1e-2) {
            grid.text("**", x, y, gp = gpar(fontsize = 8))
          } else if (facs_magma_mat[i, j] <= 0.05) {
            grid.text("*", x, y, gp = gpar(fontsize = 8))
          }
        })

##### plot correlation #####
all_results_facs = list(seismic = facs_res, scDRS = facs_scdrs_res, FUMA = facs_fuma_res, `S-MAGMA` = facs_magma_res) %>%
  map(~pivot_longer(.x, !cell_type, names_to = "trait",values_to = "Pvalue")) %>% 
  map2(names(.), ~mutate(.x, method = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>%
  group_by(trait, method) %>% 
  mutate(fdr = p.adjust(Pvalue, method="fdr")) %>%
  ungroup() %>%
  mutate(trait_name = trait_meta$official_names[match(trait,trait_meta$trait_names)]) %>%
  mutate(trait_type = trait_meta$trait_type[match(trait,trait_meta$trait_names)]) %>%
  select(cell_type,trait_name, trait_type, method, Pvalue, fdr)

cor_facs_metric = all_results_facs %>%
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

ggplot( cor_facs_metric %>%
          filter(trait_name %in% c("Bipolar disorder","Depression","Schizophrenia","Inflammatory bowel disease","Hypothroidism","Rheumatoid arthritis","Atrial fibrillation","Chronic kidney disease","LDL level","Erythrocyte count")),
        aes(x=between, y=correlation)) + 
  geom_point(size=2.5, color="#2B3990FF") +
  geom_segment(aes(xend=between, yend=0), linetype="solid",linewidth=0.75, color="#2B3990FF") +
  theme_classic() + 
  ylim(-0.1,1) + 
  coord_flip() + 
  facet_wrap(~trait_name,nrow=2) + 
  ggtitle("FACS data set, Spearman correlation across all cell types") +
  ylab("Spearman's correlation") + 
  theme(plot.title = element_text(hjust = 0.5))

ggplot( cor_facs_metric %>%
          filter(!trait_name %in% c("Bipolar disorder","Depression","Schizophrenia","Inflammatory bowel disease","Hypothroidism","Rheumatoid arthritis","Atrial fibrillation","Chronic kidney disease","LDL level","Erythrocyte count")),
        aes(x=between, y=correlation)) + 
  geom_point(size=2.5, color="#2B3990FF") +
  geom_segment(aes(xend=between, yend=0), linetype="solid",linewidth=0.75, color="#2B3990FF") +
  theme_classic() + 
  ylim(-0.1,1) + 
  coord_flip() + 
  facet_wrap(~trait_name,nrow=3) + 
  ggtitle("FACS data set, Spearman correlation across all cell types") +
  ylab("Spearman's correlation") + 
  theme(plot.title = element_text(hjust = 0.5))

##### plot difference heatmap ####
#difference heatmap function
diff_heatmap = function(fdr_mat_1, fdr_mat_2, color_vec, ...){ #neither, both, method1, method2
  ht_opt(TITLE_PADDING=unit(5,"mm"))
  agree_mat = ((fdr_mat_1<=0.05)==(fdr_mat_2<=0.05))
  color_mat = ifelse(agree_mat,
                     ifelse(fdr_mat_1<=0.05, color_vec[2],color_vec[1]),
                     ifelse(fdr_mat_1<=0.05, color_vec[3],color_vec[4]))
  #scale alpha value
  alpha_mat = ifelse(agree_mat, 1, 
                     ifelse(fdr_mat_1<=0.05, ((log(5e-2) - log(fdr_mat_1))/log(5e4))*0.9+0.1 , ((log(5e-2) - log(fdr_mat_2))/log(5e4))*0.9+0.1 ))
  for( i in 1:nrow(fdr_mat_1)){
    for(j in 1:ncol(fdr_mat_1)){
      color_mat[i,j] = scales::alpha(color_mat[i,j],alpha_mat[i, j])
    }
  }
  index_matrix = matrix(1:length(fdr_mat_1), nrow = nrow(fdr_mat_1)) %>%
    set_rownames(rownames(fdr_mat_1)) %>%
    set_colnames(colnames(fdr_mat_1))
  custom_color_function = function(index) {
    color_mat[index]
  }
  h= Heatmap(index_matrix, cluster_rows = F,cluster_columns = F,show_column_dend = F, show_row_dend = F,col=custom_color_function,show_heatmap_legend = F,
             row_names_side = "left",
             column_title_gp = gpar(fontsize=11),
             row_names_gp = gpar(fontsize=9),
             column_names_gp = gpar(fontsize=8),
             column_names_rot = 55,
             column_title_rot = 30,
             heatmap_legend_param = list(title = "association identified by"),
             rect_gp = gpar(col = "grey", lwd = 0.5),
             ...)
  color_legend = Legend(labels = names(color_vec),
                        legend_gp = gpar(fill =as.vector(color_vec)))
  col_fun_1 = colorRamp2(-log10(c(0.05,1e-6)), c(colorRampPalette(c(color_vec[3], "white"))(10)[9], color_vec[3]))
  col_fun_2 = colorRamp2(-log10(c(0.05,1e-6)), c(colorRampPalette(c(color_vec[4], "white"))(10)[9], color_vec[4]))
  alpha_legend_1 = Legend(at = -log10(c(0.05,1e-6)), col_fun = col_fun_1, labels = sprintf("%.0e", c(0.05,1e-6)), title = paste0("FDR (",names(color_vec[3]),")"))
  alpha_legend_2 = Legend(at = -log10(c(0.05,1e-6)), col_fun = col_fun_2, labels = sprintf("%.0e", c(0.05,1e-6)), title = paste0("FDR (",names(color_vec[4]),")"))
  draw(h, annotation_legend_list = list(value = color_legend, alpha = alpha_legend_1, alpha_2 = alpha_legend_2), align_heatmap_legend = "heatmap_center")
}

###scdrs with seismic
diff_heatmap(facs_fdr_mat, facs_scdrs_mat , color_vec =c( "#EEEEEE","#009E73","#D43F3AFF","#46B8DAFF" ) %>% set_names(c("neither","both","seismic","scDRS")),
             column_split =facs_ct_meta$tissue_new[match(colnames(facs_fdr_mat),facs_ct_meta$cell_type)],
             column_labels =facs_ct_meta$ontology_new[match(colnames(facs_fdr_mat),facs_ct_meta$cell_type)],
             row_split = trait_meta$trait_type[match(rownames(facs_fdr_mat), trait_meta$official_names )] %>% factor(levels = c("neuropsy","immune","others")) )

###scdrs with fuma
diff_heatmap(facs_fdr_mat, facs_fuma_mat , color_vec =c( "#EEEEEE","#009E73","#D43F3AFF","#46B8DAFF" ) %>% set_names(c("neither","both","seismic","FUMA")),
             column_split =facs_ct_meta$tissue_new[match(colnames(facs_fdr_mat),facs_ct_meta$cell_type)],
             column_labels =facs_ct_meta$ontology_new[match(colnames(facs_fdr_mat),facs_ct_meta$cell_type)],
             row_split = trait_meta$trait_type[match(rownames(facs_fdr_mat), trait_meta$official_names )] %>% factor(levels = c("neuropsy","immune","others")) )

###scdrs with magma
diff_heatmap(facs_fdr_mat, facs_magma_mat , color_vec =c( "#EEEEEE","#009E73","#D43F3AFF","#46B8DAFF" ) %>% set_names(c("neither","both","seismic","S-MAGMA")),
             column_split =facs_ct_meta$tissue_new[match(colnames(facs_fdr_mat),facs_ct_meta$cell_type)],
             column_labels =facs_ct_meta$ontology_new[match(colnames(facs_fdr_mat),facs_ct_meta$cell_type)],
             row_split = trait_meta$trait_type[match(rownames(facs_fdr_mat), trait_meta$official_names )] %>% factor(levels = c("neuropsy","immune","others")) )

##### print out tables ####
all_results_facs %>% 
  select(cell_type, trait_name, method, Pvalue, fdr) %>%
  split(.$method) %>%
  map(~select(.x, -method)) %>%
  map2(names(.), ~write.csv(.x,here("results","Tabula_muris","FACS",.y,paste0("FACS_",.y,"_association.csv")),row.names = F))
