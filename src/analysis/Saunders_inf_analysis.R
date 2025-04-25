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
library(org.Hs.eg.db)

#load data
load(here("data","expr","Saunders","Saunders_processed.rda"))

fine_cluster_res <- read.table(here("results","Saunders","fine_cluster","seismic","fine_cluster_res.txt"), header = T, sep="\t")

##### 2 calculate dfbetas #####
##calculate sscore
fine_cluster_sscore <- calc_specificity(sce = brain_sce, ct_label_col = "fine_cluster", min_avg_exp_ct = 0.01)
fine_cluster_sscore_hsa <- translate_gene_ids(fine_cluster_sscore, from = "mmu_symbol")

### for Kunkle GWAS
Kunkle_res <- fine_cluster_res %>% 
  dplyr::select(cell_type, Kunkleetal_2019) %>% 
  arrange(Kunkleetal_2019) %>% 
  mutate(FDR = p.adjust(Kunkleetal_2019, method="fdr"))

Kunkle_hc.micro_dfbetas <- find_inf_genes("HC.Microglia_Macrophage", fine_cluster_sscore_hsa, here("data","gwas","brain_gwas","zscore","Kunkleetal_2019.genes.out")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))
  
 ### for tau GWAS
tau_res <- fine_cluster_res %>% 
  dplyr::select(cell_type, tau) %>% 
  arrange(tau) %>% 
  mutate(FDR = p.adjust(tau, method="fdr"))

tau_hc.ec_dfbetas <- find_inf_genes("HC.Neurons_Entorhinal cortex", fine_cluster_sscore_hsa, here("data","gwas","brain_gwas","zscore","tau.genes.out")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))
tau_hc.ec2_dfbetas <- find_inf_genes("HC.Neurons_Medial entorhinal cortex 1", fine_cluster_sscore_hsa, here("data","gwas","brain_gwas","zscore","tau.genes.out")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))
tau_hc.ca1_dfbetas <- find_inf_genes("HC.Neurons_CA1 Principal cells", fine_cluster_sscore_hsa, here("data","gwas","brain_gwas","zscore","tau.genes.out"))  %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))

##### 3 go/kegg enrichment #####

### for Kunkle GWAS ~ HC.MG
Kunkle_hc.mg_inf_genes <- Kunkle_hc.micro_dfbetas %>%
  filter(zstat>0 & is_influential) %>%
  pull(gene)

Kunkle_hc.mg_go <- enrichGO(gene = Kunkle_hc.mg_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr",pvalueCutoff = 0.2)

Kunkle_hc.mg_kegg <- enrichKEGG(gene = Kunkle_hc.mg_inf_genes, organism = "hsa", pvalueCutoff = 0.2)

Kunkle_hc.mg_kegg.df <- Kunkle_hc.mg_kegg  %>% 
  as.data.frame() %>%
  rowwise() %>%
  mutate(new_geneID = paste0(mapIds(org.Hs.eg.db, keys= unlist(strsplit(geneID, split = "/")), keytype = "ENTREZID", column = "SYMBOL"), collapse = "/"))

### for tau GWAS ~ HC.EC / HC.EC2
tau_hc.ec_inf_genes <- tau_hc.ec_dfbetas %>%
  filter(zstat>0 & is_influential) %>%
  pull(gene)

tau_hc.ec2_inf_genes <- tau_hc.ec2_dfbetas %>%
  filter(zstat>0 & is_influential ) %>%
  pull(gene)

tau_hc.ca1_inf_genes <- tau_hc.ca1_dfbetas %>%
  filter(zstat>0 & is_influential ) %>%
  pull(gene)

tau_hc.ec_go <- enrichGO(gene = tau_hc.ec_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr", pvalueCutoff = 0.2)
tau_hc.ec2_go <- enrichGO(gene = tau_hc.ec2_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr", pvalueCutoff = 0.2)
tau_hc.ca1_go <- enrichGO(gene = tau_hc.ca1_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr", pvalueCutoff = 0.2)


tau_hc.ec_kegg <- enrichKEGG(gene =  tau_hc.ec_inf_genes, organism = "hsa", pvalueCutoff = 0.2)
tau_hc.ec2_kegg <- enrichKEGG(gene =  tau_hc.ec2_inf_genes, organism = "hsa", pvalueCutoff = 0.2)
tau_hc.ca1_kegg <- enrichKEGG(gene =  tau_hc.ca1_inf_genes, organism = "hsa", pvalueCutoff = 0.2)

tau_hc.ec_kegg.df <- tau_hc.ec_kegg %>% 
  as.data.frame() %>%
  rowwise() %>%
  mutate(new_geneID = paste0(mapIds(org.Hs.eg.db, keys= unlist(strsplit(geneID, split = "/")), keytype = "ENTREZID", column = "SYMBOL"), collapse = "/"))

tau_hc.ec2_kegg.df <- tau_hc.ec2_kegg %>% 
  as.data.frame() %>%
  rowwise() %>%
  mutate(new_geneID = paste0(mapIds(org.Hs.eg.db, keys= unlist(strsplit(geneID, split = "/")), keytype = "ENTREZID", column = "SYMBOL"), collapse = "/"))

tau_hc.ca1_kegg.df <- tau_hc.ca1_kegg %>% 
  as.data.frame() %>%
  rowwise() %>%
  mutate(new_geneID = paste0(mapIds(org.Hs.eg.db, keys= unlist(strsplit(geneID, split = "/")), keytype = "ENTREZID", column = "SYMBOL"), collapse = "/"))

##### 4 save results ####
#mircoglia ~ clinical AD
Kunkle_hc.micro_dfbetas %>%
  relocate(gene, genename) %>%
  set_colnames(c("entrez id", "gene symbol", "specificity score", "gene risk (z-score)", "dfbetas", "influential")) %>%
  write.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_dfbetas.txt"), row.names = F, quote = F, sep="\t")

Kunkle_hc.mg_go %>% 
  as.data.frame() %>% 
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "GO name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_go.txt"), row.names = F, quote = F,sep="\t")

Kunkle_hc.mg_kegg.df %>% 
  ungroup() %>%
  dplyr::select(ID, Description, p.adjust, Count, new_geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "KEGG name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_kegg.txt"),row.names = F, quote = F, sep="\t")

#ec ~ tau
tau_hc.ec_dfbetas %>% 
  relocate(gene, genename) %>%
  set_colnames(c("entrez id", "gene symbol", "specificity score", "gene risk (z-score)", "dfbetas", "influential")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ec_dfbetas.txt"),row.names = F, quote = F, sep="\t")

tau_hc.ec_go %>% 
  as.data.frame() %>% 
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "GO name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ec_go.txt"),row.names = F, quote = F, sep="\t")

tau_hc.ec_kegg.df %>%  
  ungroup() %>%
  dplyr::select(ID, Description, p.adjust, Count, new_geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "KEGG name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ec_kegg.txt"),row.names = F, quote = F, sep="\t")

#ec2 ~ tau
tau_hc.ec2_dfbetas %>% 
  relocate(gene, genename) %>%
  set_colnames(c("entrez id", "gene symbol", "specificity score", "gene risk (z-score)", "dfbetas", "influential")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ec2_dfbetas.txt"),row.names = F, quote = F,  sep="\t")

tau_hc.ec2_go %>% 
  as.data.frame() %>% 
  filter(Count >= 3) %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  set_colnames(c("ID", "GO name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ec2_go.txt"),row.names = F, quote = F, sep="\t")

tau_hc.ec2_kegg.df %>%  
  ungroup() %>%
  dplyr::select(ID, Description, p.adjust,Count, new_geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "KEGG name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ec2_kegg.txt"),row.names = F, quote = F, sep="\t")

#ca1 ~ tau
tau_hc.ca1_dfbetas %>% 
  relocate(gene, genename) %>%
  set_colnames(c("entrez id", "gene symbol", "specificity score", "gene risk (z-score)", "dfbetas", "influential")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ca1_dfbetas.txt"), row.names = F, quote = F, sep="\t")

tau_hc.ca1_go %>% 
  as.data.frame() %>% 
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "GO name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ca1_go.txt"), row.names = F, quote = F, sep="\t")

tau_hc.ca1_kegg.df %>%  
  ungroup() %>%
  dplyr::select(ID, Description, p.adjust,Count, new_geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "KEGG name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results","Saunders","inf_analysis","tau_hc_ca1_kegg.txt"), row.names = F, quote = F, sep="\t")

