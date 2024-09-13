#plot for the simulation results
##### load packages and results ####
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

###load previously stored metadata
#only keep cell types with single mapping to cell ontology and tissue
facs_cl_mapping <- read.table(here("results","Tabula_muris","FACS","facs_ontology.txt"),sep="\t", header = T) %>% 
  as_tibble() %>%
  group_by(cell_ontology_id, tissue) %>%
  filter(n()<2) %>%
  ungroup
droplet_cl_mapping <- read.table(here("results","Tabula_muris","droplet","droplet_ontology.txt"),sep="\t", header = T) %>% 
  as_tibble() %>%
  group_by(cell_ontology_id, tissue) %>%
  filter(n()<2) %>%
  ungroup
ts_cl_mapping <- read.table(here("results","Tabula_sapiens","ts_ontology.txt"),sep="\t", header = T) %>% 
  as_tibble() %>%
  group_by(cell_ontology_id, tissue) %>%
  filter(n()<2) %>%
  ungroup

###load results
facs_res <- read.table(here("results","Tabula_muris","FACS","seismic","new_facs_res.txt"), sep="\t", header = T) %>% 
  as_tibble() %>%
  pivot_longer(!cell_type, names_to = "trait", values_to = "Pvalue")
droplet_res <- read.table(here("results","Tabula_muris","droplet","seismic","new_droplet_res.txt"), sep="\t", header = T) %>% 
  as_tibble() %>%
  pivot_longer(!cell_type, names_to = "trait", values_to = "Pvalue")
ts_res <- read.table(here("results","Tabula_sapiens","seismic","new_ts_res.txt"), sep="\t", header = T) %>% 
  as_tibble() %>%
  pivot_longer(!cell_type, names_to = "trait", values_to = "Pvalue")

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


##### extract common pairs ####
facs_droplet_common <- facs_cl_mapping %>%  
  inner_join(droplet_cl_mapping %>% mutate(tissue = recode(tissue, "Heart_and_Aorta" = "Heart")) ,  #recode tissue names for mapping
            by = c("cell_ontology_id" = "cell_ontology_id" , "tissue" = "tissue")) %>%
  mutate(unified_id = paste0(tissue,".",cell_ontology_id))

facs_ts_common <- facs_cl_mapping %>%  
  inner_join(ts_cl_mapping %>% mutate(tissue = recode(tissue,"Mammary" = "Mammary_Gland", "Muscle" = "Limb_Muscle", "Bone_Marrow" = "Marrow")) ,
            by = c("cell_ontology_id" = "cell_ontology_id" , "tissue" = "tissue")) %>%
  mutate(unified_id = paste0(tissue,".",cell_ontology_id))

facs_droplet_pairs <- facs_droplet_common %>% 
  left_join(facs_res, by = c("cluster_name.x" = "cell_type")) %>%
  left_join(droplet_res, by = c("cluster_name.y" = "cell_type", "trait"="trait")) %>% 
  drop_na() %>% #drop ones without appearance in facs/droplet
  select(unified_id, trait, Pvalue.x, Pvalue.y)
  
facs_ts_pairs <- facs_ts_common %>% 
  left_join(facs_res, by = c("cluster_name.x" = "cell_type")) %>%
  left_join(ts_res, by = c("cluster_name.y" = "cell_type", "trait"="trait")) %>% 
  select(unified_id, trait, Pvalue.x, Pvalue.y) %>%
  drop_na() #drop ones without appearance in facs/dts

##### calculate correlation ###### and make the plot
facs_droplet_cor_df <- facs_droplet_pairs %>% 
  mutate_if(is.numeric, ~-log10(.)) %>%
  group_by(trait) %>%
  mutate(correlation = cor(Pvalue.x, Pvalue.y ,method="spearman")) %>% 
  ungroup() %>%
  distinct(trait, correlation) %>%
  left_join(trait_meta, by = c("trait"="trait_names")) %>%
  dplyr::rename("trait_names" = "official_names") %>%
  mutate(trait_names = factor(trait_names, levels=rev(trait_meta$official_names)))

facs_ts_cor_df <- facs_ts_pairs %>% 
  mutate_if(is.numeric, ~-log10(.)) %>%
  group_by(trait) %>%
  mutate(correlation = cor(Pvalue.x, Pvalue.y ,method="spearman")) %>% 
  ungroup() %>%
  distinct(trait, correlation) %>%
  left_join(trait_meta, by = c("trait"="trait_names")) %>%
  dplyr::rename("trait_names" = "official_names") %>%
  mutate(trait_names = factor(trait_names, levels=rev(trait_meta$official_names)))



##plot
ggplot( rbind(facs_droplet_cor_df %>% mutate(between = "TM FACS and TM droplet"), facs_ts_cor_df %>% mutate(between = "TM FACS and TS")), aes(x=trait_names, y=correlation,color=trait_type)) + 
  geom_point(size=2.5,alpha=0.8) +
  geom_segment(aes(xend=trait_names, yend=0), linetype="solid",linewidth=0.75,alpha=0.8) +
  theme_classic() + 
  ylim(0,1) + 
  coord_flip() + 
  facet_wrap(~between,nrow=1,scales="free_y") +
  ggtitle("Correlation of cell type association P-values between datasets") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#566371FF","#C6DC65FF","#C54D58FF"))
