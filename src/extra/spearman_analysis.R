#Analysis of Tabula muris data set using Spearman's correlation model
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

##### calculate specificity score and enrichment results #####
facs_sscore  <- calc_specificity(sce = facs_obj_sce , ct_label_col = "cluster_name")
droplet_sscore  <- calc_specificity(sce = droplet_obj_sce , ct_label_col = "cluster_name", min_avg_exp_ct = 0.01) #for droplet

#map to human genes
facs_sscore_hsa <- translate_gene_ids(facs_sscore, from = "mmu_symbol")
droplet_sscore_hsa <- translate_gene_ids(droplet_sscore, from = "mmu_symbol")

##enrichment
#magma zscore file
magma_zscore_file <- list.files(here("data","gwas","tm_gwas","zscore"), full.names = T) %>%
  set_names(str_extract(string = ., pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.genes\\.out", replacement = "", x=.))

#enrichment
facs_association_spearman <- magma_zscore_file %>%
  map(~get_ct_trait_associations(sscore = facs_sscore_hsa, magma = .x, model = "spearman"))
droplet_association_spearman <- magma_zscore_file %>%
  map(~get_ct_trait_associations(sscore = droplet_sscore_hsa, magma = .x, model = "spearman"))

##save results
facs_res_spearman <- facs_association_spearman %>%
  map(~select(.x, cell_type, pvalue)) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by="cell_type"))

droplet_res_spearman <- droplet_association_spearman %>%
  map(~select(.x, cell_type, pvalue)) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by="cell_type"))

#
write.table(facs_res_spearman, here("results","Tabula_muris","FACS","seismic","facs_res.spearman.txt"),quote=F, sep="\t", row.names = F)
write.table(droplet_res_spearman, here("results","Tabula_muris","droplet","seismic","droplet_res.spearman.txt"),quote=F, sep="\t", row.names = F)


##### Load previous results #####
neuropsy_disesaes <- c("insominia","BD","MDD","Nrt_new","Scz_new")
immune_diseases <- c("AID","Hypothyroidism","RA","ibd","cd","uc")
others <- c("Smoking","BMI","College_edu","RBC","Lymphocyte_count","Monocyte_count","T1D","T2D_2","Cardiovas","SBP","AF","CKD","glucose_q","HDL_q","LDL_q","TG_q")
trait_meta <- tibble(trait_names = c(neuropsy_disesaes, immune_diseases, others),
                     official_names = c("Insomnia","Bipolar disorder","Depression","Neuroticism","Schizophrenia",
                                        "Autoimmune diseases","Hypothroidism","Rheumatoid arthritis","Inflammatory bowel disease","Crohn's disease","Ulcerative colitis",
                                        "Smoking","BMI","College education","Erythrocyte count","Lymphocyte count","Monocyte count","Type I diabetes","Type II diabetes","Cardiovascular diseases",
                                        "Systolic blood pressure","Atrial fibrillation","Chronic kidney disease","Glucose level","HDL level","LDL level","Triglyceride level"),
                     trait_type = factor(c(rep("neuropsy",length(neuropsy_disesaes)), rep("immune",length(immune_diseases)), rep("others",length(others))), levels=c("neuropsy","immune","others")) )

#load results
facs_res <- read.table( here("results","Tabula_muris","FACS","seismic","facs_res.txt"), header = T, sep = "\t") %>% 
  as_tibble() %>%
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))  %>%
  mutate(cell_type = factor(cell_type, levels = sort(cell_type))) %>%
  arrange(cell_type)

droplet_res <- read.table( here("results","Tabula_muris","droplet","seismic","droplet_res.txt"), header = T, sep = "\t") %>% 
  as_tibble() %>%
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))  %>%
  mutate(cell_type = factor(cell_type, levels = sort(cell_type))) %>%
  arrange(cell_type)

#calculate correlation
facs_all <- facs_res_spearman %>% 
  pivot_longer(cols = !cell_type, names_to = "trait", values_to = "facs_spearman") %>%
  left_join(pivot_longer(facs_res, cols = !cell_type, names_to = "trait", values_to = "facs_linear"), by = c("cell_type","trait")) 
  
droplet_all <- droplet_res_spearman %>% 
  pivot_longer(cols = !cell_type, names_to = "trait", values_to = "droplet_spearman") %>%
  left_join(pivot_longer(droplet_res, cols = !cell_type, names_to = "trait", values_to = "droplet_linear"), by = c("cell_type","trait")) 

facs_cor <- facs_all %>%
  group_by(trait) %>%
  summarise(correlation = cor(-log(facs_spearman), -log(facs_linear))) %>%
  mutate(data_set = "TM FACS")

droplet_cor <- droplet_all %>%
  group_by(trait) %>%
  summarise(correlation = cor(-log(droplet_spearman), -log(droplet_linear))) %>%
  mutate(data_set = "TM droplet")

all_data <- rbind(facs_cor, droplet_cor) 

#plot
all_data <- all_data %>%
  left_join(trait_meta, by = c("trait" = "trait_names")) %>%
  mutate(official_names = factor(official_names, levels = trait_meta$official_names))

ggplot(all_data, aes(x=official_names, y=correlation,color=trait_type)) + 
  geom_point(size=2.5,alpha=0.8) +
  geom_segment(aes(xend=official_names, yend=0), linetype="solid",linewidth=0.75,alpha=0.8) +
  theme_classic() + 
  facet_wrap(~data_set) +
  ylim(0,1) + 
  coord_flip() + 
  ggtitle("Mean correlation of cell type association P-values between Sperman's correlation model and standard linear model") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#566371FF","#C6DC65FF","#C54D58FF"))
