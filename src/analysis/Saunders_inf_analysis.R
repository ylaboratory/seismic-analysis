#Analysis of Saunders et al data set
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

if (!require("anndata")){
  install.packages("anndata")
  library("anndata")
}

if (!require("clusterProfiler")){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("clusterProfiler")
  library("clusterProfiler")
} 
if (!require("org.Hs.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("org.Hs.eg.db")
  library("org.Hs.eg.db")
} 

##### 2 calculate dfbetas #####
### for Kunkle GWAS
Kunkle_res = fine_cluster_res %>% 
  select(cell_type, Kunkleetal_2019) %>% 
  arrange(Kunkleetal_2019) %>% 
  mutate(FDR = p.adjust(Kunkleetal_2019, method="fdr"))

Kunkle_hc.micro_dfbetas = gene_inf_measure(brain_sce, gene_zscore_df = gwas_zscore_AD, cell_type = "HC.Microglia_Macrophage",trait_name = "Kunkleetal_2019") 

### for tau GWAS
tau_res = fine_cluster_res %>% 
  select(cell_type, tau) %>% 
  arrange(tau) %>% 
  mutate(FDR = p.adjust(tau, method="fdr"))

tau_hc.ec_dfbetas = gene_inf_measure(brain_sce, gene_zscore_df = gwas_zscore_AD, cell_type = "HC.Neurons_Entorhinal cortex",trait_name = "tau") 
tau_hc.ec2_dfbetas = gene_inf_measure(brain_sce, gene_zscore_df = gwas_zscore_AD, cell_type = "HC.Neurons_Medial entorhinal cortex 1",trait_name = "tau") 

##### 3 go enrichment #####
### for Kunkle GWAS ~ HC.MG
Kunkle_hc.mg_inf_genes = Kunkle_hc_micro_dfbetas %>%
  filter(Kunkleetal_2019_z_stat>0 & influential) %>%
  pull(hsa_entrez)

Kunkle_hc.mg_go = enrichGO(gene = Kunkle_HC.mg_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")

### for tau GWAS ~ HC.EC / HC.EC2
tau_hc.ec_inf_genes = tau_hc.ec_dfbetas %>%
  filter(tau_z_stat>0 & influential) %>%
  pull(hsa_entrez)
tau_hc.ec2_inf_genes = tau_hc.ec2_dfbetas %>%
  filter(tau_z_stat>0 & influential) %>%
  pull(hsa_entrez)

tau_hc.ec_go = enrichGO(gene = tau_HC.ec_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")
tau_hc.ec2_go = enrichGO(gene = tau_HC.ec2_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")

##### 4 save results ####
Kunkle_hc_micro_dfbetas %>% write.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_dfbetas.txt"),row.names = F,quote=F, sep="\t")
Kunkle_hc.mg_go %>% as.data.frame() %>% write.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_go.txt"),row.names = F,quote=F, sep="\t")

tau_hc.ec_dfbetas %>% write.table(here("results","Saunders","inf_analysis","Kunkle_hc_ec_dfbetas.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ec2_dfbetas %>% write.table(here("results","Saunders","inf_analysis","Kunkle_hc_ec2_dfbetas.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ec_go %>% as.data.frame() %>% write.table(here("results","Saunders","inf_analysis","Kunkle_hc_ec_go.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ec2_go %>% as.data.frame() %>% write.table(here("results","Saunders","inf_analysis","Kunkle_hc_ec2_go.txt"),row.names = F,quote=F, sep="\t")