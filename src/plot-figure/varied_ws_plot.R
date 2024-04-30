#generate plots for results of the varied window size plot
##### load packages and results #####
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

#load trait meta information
neuropsy_disesaes <- c("insominia","BD","MDD","Nrt_new","Scz_new")
immune_diseases <- c("AID","Hypothyroidism","RA","ibd","cd","uc")
others <- c("Smoking","BMI","College_edu","RBC","Lymphocyte_count","Monocyte_count","T1D","T2D_2","Cardiovas","SBP","AF","CKD","glucose_q","HDL_q","LDL_q","TG_q")
trait_meta <- tibble(trait_names = c(neuropsy_disesaes, immune_diseases, others),
                    official_names = c("Insomnia","Bipolar disorder","Depression","Neuroticism","Schizophrenia",
                                       "Autoimmune diseases","Hypothroidism","Rheumatoid arthritis","Inflammatory bowel disease","Crohn's disease","Ulcerative colitis",
                                       "Smoking","BMI","College education","Erythrocyte count","Lymphocyte count","Monocyte count","Type I diabetes","Type II diabetes","Cardiovascular diseases",
                                       "Systolic blood pressure","Atrial fibrillation","Chronic kidney disease","Glucose level","HDL level","LDL level","Triglyceride level"),
                    trait_type = factor(c(rep("neuropsy",length(neuropsy_disesaes)), rep("immune",length(immune_diseases)), rep("others",length(others))), levels=c("neuropsy","immune","others")) )

varied_ws_res <- list.files(here("results","varied_ws"), pattern="new_all_res") %>%
  set_names(gsub(pattern = "new_all_res\\.|\\.txt",  x= unlist(.), replacement = "")) %>%
  map(~read.table(here("results","varied_ws",.x), header = T, sep = "\t")) %>% 
  map(~pivot_longer(.x, !cell_type, names_to = "trait",values_to = "Pvalue")) %>%
  map2(names(.), ~mutate(.x, window_range = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>%
  mutate(Pvalue = -log10(Pvalue)) #-log10 Pvalues

varied_ws_cor <- varied_ws_res %>% 
  split(.$trait) %>% 
  map(~select(.x, - trait)) %>% 
  map(~pivot_wider(.x, names_from = window_range, values_from = Pvalue)) %>%
  map(~select(.x, - cell_type)) %>%
  map(~cor(.x)) %>% #correlation
  map(~reshape2::melt(.x)) %>%  #remodel shape
  map2(names(.), ~mutate(.x, trait = .y)) %>% #combing with trait names
  purrr::reduce(~rbind(.x,.y)) %>%
  mutate(trait=trait_meta$official_names[match(trait,trait_meta$trait_names)]) %>%
  mutate(trait=factor(trait,levels=trait_meta$official_names)) %>%
  mutate(Var1 = factor(Var1, levels = c("10_1.5","10_10","20_5","35_10","50_10","50_50","100_10","100_50")),
         Var2 = factor(Var2, levels = c("10_1.5","10_10","20_5","35_10","50_10","50_50","100_10","100_50")))


###### plot results #####                                 
ggplot(varied_ws_cor, aes(Var1,Var2,fill=value)) + 
  geom_tile(color="grey60")+
  facet_wrap(~ trait) +
  theme_minimal() +
  labs(x = 'MAGMA gene analysis window size: downstream(kb)_upstream(kb)', y = 'MAGMA gene analysis window size: downstream(kb)_upstream(kb)',fill="correlation") +
  ggtitle("Correlation of results by semsmicGWAS with varying window sizes")+
  viridis::scale_fill_viridis() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1),plot.title = element_text(hjust = 0.5))
