#Analysis of Tabula muris data set - split the data by donor ID
##### 1. load packages and data#######
###load packages
if (!require("here")) {
  install.packages("here")
  library("here")
}

if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}

if (!require("seismicGWAS")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

###load data
load(here("data","expr","Tabula_muris","TM_processed.rda"))

##### 1. separate cells by individuals #####
facs_obj_ind <- unique(facs_obj_sce$mouse.id) %>%
  set_names(.) %>%
  map(~facs_obj_sce[, facs_obj_sce$mouse.id == .x]) %>%
  keep(~ncol(.) > 1000)


##### 2. calculate specificity score and associations #####
facs_sscore_ind <- facs_obj_ind %>%
  map(~calc_specificity(sce = .x , ct_label_col = "cluster_name")) %>% #calculate specificity score
  map(~translate_gene_ids(.x, from = "mmu_symbol"))  #map to hsa

#get the number of cell types
max_num_ct <- facs_sscore_ind %>% 
  map(~ncol(.x)) %>%
  unlist() %>%
  max

facs_sscore_ind <- facs_sscore_ind %>% 
  keep(~ncol(.x) > max_num_ct/2) %>%
  map(~as.matrix(.x))


#magma zscore file
magma_zscore_file <- list.files(here("data","gwas","tm_gwas","zscore"), full.names = T) %>%
  set_names(str_extract(string = ., pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.genes\\.out", replacement = "", x=.))

#get associations
facs_ind_association <- magma_zscore_file %>%
  map(~ {zscore_file <- .x; map(facs_sscore_ind, ~get_ct_trait_associations(.x, magma = zscore_file))}) 


##### 2.5 Calculate the correlation of specificity score #####
all_sscore_5 <- facs_obj_ind[names(facs_sscore_ind)] %>%
  map(~calc_specificity(sce = .x , ct_label_col = "cluster_name",min_ct_size = 5)) %>% #calculate specificity score
  map(~translate_gene_ids(.x, from = "mmu_symbol"))

facs_score_cor <- all_sscore_5 %>%
  map(~colnames(.x)) %>% 
  purrr::reduce(~c(.x, .y)) %>%
  unique() %>%
  set_names(.) %>%
  map(~ {ct = .x; map(all_sscore_5, ~.x[,which(colnames(.x) == ct)])}) %>% 
  map(~map(.x, ~tibble(gene = names(.x), value = as.vector(.x)))) %>%
  map(~map2(.x, names(.x), ~set_colnames(.x, c("gene", .y)))) %>%
  map(~keep(.x, ~ncol(.) > 1)) %>%
  map(~purrr::reduce(.x, ~left_join(.x, .y, by = "gene"))) %>%
  map(~select(.x, -gene)) %>%
  keep(~ncol(.) >= 2) %>%
  map(~cor(.x, use = "na.or.complete")) %>% 
  map(~(sum(.x) - nrow(.x))/( nrow(.x)^2 - nrow(.x)))

num_ct <- map(names(facs_score_cor), ~length(which(facs_obj_sce$cluster_name == .x)))

num_ind <- map(names(facs_score_cor), ~facs_obj_sce$mouse.id[facs_obj_sce$cluster_name == .x]) %>%
  map(~ {tot_id = .x; map(names(all_sscore_5), ~length(which(tot_id == .x)))}) %>%
  map(~keep(.x, .>=3)) %>%
  map(~min(unlist(.x)))

score_cor_df <- tibble(score_correlation = unlist(facs_score_cor), 
                       cell_type = names(facs_score_cor), 
                       tot_cell_num = unlist(num_ct), 
                       min_cell_number_ind = unlist(num_ind),
                       filtered_out = ifelse(tot_cell_num < 20, T, F))

ggplot(score_cor_df, aes(x = min_cell_number_ind, y = score_correlation, color = log10(tot_cell_num), shape = filtered_out)) + 
  geom_point(size=3, alpha = 0.8) +
  scale_color_viridis_b() +
  theme_classic() +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_log10() +
  #geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("minimum number of cells observed in an donor") +
  ylab("mean correlation across donor pairs")  +
  guides(fill = guide_legend(title="log10(tot number of cells)")) 

##### 3. Calculate the correlation for results#####
#reshape format
facs_ind_association <- facs_ind_association %>% 
  map(~map(.x, ~select(.x, -FDR))) %>%
  map(~map2(.x,names(.x), ~set_colnames(.x, c("cell_type", .y)))) %>%
  map(~purrr::reduce(.x, ~full_join(.x, .y, by = "cell_type")))  


#calculate the correlation matrix
facs_ind_correlation <- facs_ind_association %>%
  map(~select(.x, - cell_type)) %>%
  map(~mutate_if(.x, is.numeric, ~-log10(.))) %>% 
  map(~cor(.x, use = "na.or.complete"))



##### 4. plot ####
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

facs_cor_df <- facs_ind_correlation %>%
  map(~(sum(.x) - nrow(.x))/(nrow(.x)*(nrow(.x)-1))) %>%
  unlist() %>%
  as.data.frame() %>% 
  as_tibble(rownames = "trait") %>%
  set_colnames(c("trait","correlation")) %>%
  left_join(trait_meta, by = c("trait" = "trait_names")) %>%
  mutate(official_names = factor(official_names, levels = rev(trait_meta$official_names))) %>%
  arrange(official_names)

ggplot(facs_cor_df, aes(x=official_names, y=correlation,color=trait_type)) + 
  geom_point(size=2.5,alpha=0.8) +
  geom_segment(aes(xend=official_names, yend=0), linetype="solid",linewidth=0.75,alpha=0.8) +
  theme_classic() + 
  ylim(0,1) + 
  coord_flip() + 
  ggtitle("Mean correlation of cell type association P-values between all indvidual pairs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#566371FF","#C6DC65FF","#C54D58FF"))
