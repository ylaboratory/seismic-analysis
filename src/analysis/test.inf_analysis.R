#Analysis of Tabula muris data set
##### 1. load packages and data#######
###load packages
if (!require("here")){
  install.packages("here")
  library("here")
}

if (!require("magrittr")){
  install.packages("magrittr")
  library("magrittr")
}
if (!require("tidyverse")){
  install.packages("tidyverse")
  library("tidyverse")
}

if (!require("seismicGWAS")){
  if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}
library(org.Hs.eg.db)
library(clusterProfiler)
#import results
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

#load results - facs
facs_res <- read.table( here("results","Tabula_muris","FACS","seismic","new_facs_res.txt"), header = T, sep = "\t") %>% 
  as_tibble() %>%
  select(all_of(c("cell_type", neuropsy_disesaes, immune_diseases ,others)))  %>%
  mutate(cell_type = factor(cell_type, levels = sort(cell_type))) %>%
  arrange(cell_type)


##sscore
load(here("data","expr","Tabula_muris","TM_processed.rda"))

facs_sscore  <- calc_specificity(sce = facs_obj_sce , ct_label_col = "cluster_name")

facs_sscore_hsa <- translate_gene_ids(facs_sscore, from = "mmu_symbol")

##influential gene analysis
#for SBP
SBP.endo_dfbetas <- find_inf_genes("Heart.endothelial cell", facs_sscore_hsa, here("data","gwas","tm_gwas","zscore","SBP.genes.out"))
SBP.smc_dfbetas <- find_inf_genes("Heart.smooth muscle cell", facs_sscore_hsa, here("data","gwas","tm_gwas","zscore","SBP.genes.out"))
SBP.mf_dfbetas <- find_inf_genes("Heart.myofibroblast cell", facs_sscore_hsa, here("data","gwas","tm_gwas","zscore","SBP.genes.out"))

SBP.endo_genes <- SBP.endo_dfbetas %>% filter(is_influential) %>% pull(gene)
SBP.smc_genes <- SBP.smc_dfbetas %>% filter(is_influential) %>% pull(gene)
SBP.mf_genes <- SBP.mf_dfbetas %>% filter(is_influential) %>% pull(gene)

SBP.endo_go <- enrichGO(gene =SBP.endo_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")
SBP.smc_go <- enrichGO(gene =SBP.smc_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")
SBP.mf_go <- enrichGO(gene =SBP.mf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")

#for LDL
LDL.hpt_dfbetas <- find_inf_genes("Liver.hepatocyte", facs_sscore_hsa, here("data","gwas","tm_gwas","zscore","LDL_q.genes.out"))
LDL.ent_dfbetas <- find_inf_genes("Large_Intestine.Enterocyte (Proximal)", facs_sscore_hsa, here("data","gwas","tm_gwas","zscore","LDL_q.genes.out"))
LDL.b_dfbetas <- find_inf_genes("Liver.B cell", facs_sscore_hsa, here("data","gwas","tm_gwas","zscore","LDL_q.genes.out"))

LDL.hpt_genes <- LDL.hpt_dfbetas %>% filter(is_influential) %>% pull(gene)
LDL.ent_genes <- LDL.ent_dfbetas %>% filter(is_influential) %>% pull(gene)
LDL.b_genes <- LDL.b_dfbetas %>% filter(is_influential) %>% pull(gene)

LDL.hpt_go <- enrichGO(gene =LDL.hpt_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")
LDL.ent_go <- enrichGO(gene =LDL.ent_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")
LDL.b_go <- enrichGO(gene =LDL.ent_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")

#for 
