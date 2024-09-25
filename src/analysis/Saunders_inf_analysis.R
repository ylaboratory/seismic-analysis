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
  select(cell_type, tau) %>% 
  arrange(tau) %>% 
  mutate(FDR = p.adjust(tau, method="fdr"))

tau_hc.ec_dfbetas <- find_inf_genes("HC.Neurons_Entorhinal cortex", fine_cluster_sscore_hsa, here("data","gwas","brain_gwas","zscore","tau.genes.out")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))
tau_hc.ec2_dfbetas <- find_inf_genes("HC.Neurons_Medial entorhinal cortex 1", fine_cluster_sscore_hsa, here("data","gwas","brain_gwas","zscore","tau.genes.out")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))
tau_hc.ca1_dfbetas <- find_inf_genes("HC.Neurons_CA1 Principal cells", fine_cluster_sscore_hsa, here("data","gwas","brain_gwas","zscore","tau.genes.out"))  %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))

##### 3 go enrichment #####

### for Kunkle GWAS ~ HC.MG
Kunkle_hc.mg_inf_genes <- Kunkle_hc.micro_dfbetas %>%
  filter(zstat>0 & is_influential) %>%
  pull(gene)

Kunkle_hc.mg_go <- enrichGO(gene = Kunkle_hc.mg_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr",pvalueCutoff = 0.2)
Kunkle_hc.mg_kegg <- enrichKEGG(gene = Kunkle_hc.mg_inf_genes, organism = "hsa", pvalueCutoff = 0.2)

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

#t = read.table("results/Saunders/inf_analysis/Kunkle_hc_ec_dfbetas.txt", header=T)
#t2 = t %>% filter(tau_z_stat>0, influential) %>% pull(hsa_entrez)
#t3 = read.table("results/Saunders/inf_analysis/Kunkle_hc_ec2_dfbetas.txt", header=T)
#t4 = t3 %>% filter(tau_z_stat>0, influential)

#tau_hc.ec_go <- enrichGO(gene =  t2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")
#tau_hc.ec2_go <- enrichGO(gene = t4 , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr")

tau_hc.ec_kegg <- enrichKEGG(gene =  tau_hc.ec_inf_genes, organism = "hsa", pvalueCutoff = 0.2)
tau_hc.ec2_kegg <- enrichKEGG(gene =  tau_hc.ec2_inf_genes, organism = "hsa", pvalueCutoff = 0.2)
tau_hc.ca1_kegg <- enrichKEGG(gene =  tau_hc.ca1_inf_genes, organism = "hsa", pvalueCutoff = 0.2)

tau_hc.overlap_kegg <- enrichKEGG(gene = intersect(tau_hc.ec_inf_genes, tau_hc.ec2_inf_genes), organism = "hsa",pvalueCutoff = 0.2)

###plot 
### kegg plot
jaccard_idx <- function(geneid1, geneid2){
  gene_id1 <- strsplit(geneid1, split="/") %>% unlist
  gene_id2 <- strsplit(geneid2, split="/") %>% unlist
  return(length(intersect(gene_id1, gene_id2))/length(union(gene_id1, gene_id2)))
}
Kunkle_hc.mg_kegg.df <- as.data.frame(Kunkle_hc.mg_kegg) %>%
  rowwise() %>%
  mutate(geneNames = paste0(mapIds(org.Hs.eg.db, keys = unlist(strsplit(geneID, split="/")), column = "SYMBOL", keytype = "ENTREZID"), collapse = "/"))

tau_hc.ec_kegg.df <- as.data.frame(tau_hc.ec_kegg) %>%
  rowwise() %>%
  mutate(geneNames = paste0(mapIds(org.Hs.eg.db, keys = unlist(strsplit(geneID, split="/")), column = "SYMBOL", keytype = "ENTREZID"), collapse = "/"))

tau_hc.ec2_kegg.df <- as.data.frame(tau_hc.ec2_kegg) %>%
  rowwise() %>%
  mutate(geneNames = paste0(mapIds(org.Hs.eg.db, keys = unlist(strsplit(geneID, split="/")), column = "SYMBOL", keytype = "ENTREZID"), collapse = "/"))
  
tau_hc.ca1_kegg.df <- as.data.frame(tau_hc.ca1_kegg) %>%
  rowwise() %>%
  mutate(geneNames = paste0(mapIds(org.Hs.eg.db, keys = unlist(strsplit(geneID, split="/")), column = "SYMBOL", keytype = "ENTREZID"), collapse = "/"))

#plot venn diagram
library(ggvenn)
gene_venn <- list(EC = tau_hc.ec_inf_genes , EC2 = tau_hc.ec2_inf_genes, CA1 = tau_hc.ca1_inf_genes, MG = Kunkle_hc.mg_inf_genes) 
ggvenn(gene_venn,text_size = 6,set_name_size=6,show_percentage = F,fill_alpha = 0.4) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-2, 1.5)) 

kegg_venn <- list(EC = tau_hc.ec_kegg.df$ID , EC2 = tau_hc.ec2_kegg.df$ID , CA1 = NULL, MG = Kunkle_hc.mg_kegg.df$ID) 
ggvenn(kegg_venn,text_size = 6,set_name_size=6,show_percentage = F,fill_alpha = 0.4) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-2, 1.5))

#plot a network 
kegg_data <- rbind(tau_hc.ec_kegg.df, tau_hc.ec2_kegg.df) %>%
  as_tibble() %>%
  mutate(cell_type = c(rep("EC", nrow(tau_hc.ec_kegg)), rep("EC2", nrow(tau_hc.ec2_kegg)) )) 
kegg_data_nw <- map(kegg_data$geneID, 
                     ~{
                       gene_id1 <- .x
                       map(kegg_data$geneID, ~jaccard_idx(gene_id1, .x))
                     }) %>%
  unlist() %>%
  matrix(nrow = nrow(kegg_data))

library(ComplexHeatmap)
library(circlize)
ht_opt(TITLE_PADDING=unit(5,"mm"))
col_fun <- colorRamp2(c(0, 0.33, 0.66,1), viridis::viridis(4))

kegg_data_nw <- kegg_data_nw %>% 
  set_rownames( kegg_data$Description) %>%
  set_colnames( kegg_data$Description) 

Heatmap(kegg_data_nw, name="Jaccard idx",cluster_rows = T,cluster_columns = T,show_column_dend = F, show_row_dend = F, rect_gp = gpar(col = "grey", lwd = 0.5), col = col_fun, 
        column_split = kegg_data$cell_type, 
        row_split = kegg_data$cell_type,
        row_names_side = "left",
        column_title_gp = gpar(fontsize=11),
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=7),
        column_names_rot = 55,
        column_title_rot = 30)

##### 4 save results ####
Kunkle_hc.micro_dfbetas %>% write.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_dfbetas.txt"),row.names = F,quote=F, sep="\t")
Kunkle_hc.mg_go %>% as.data.frame() %>% write.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_go.txt"),row.names = F,quote=F, sep="\t")
Kunkle_hc.mg_kegg.df %>% 
  ungroup() %>%
  mutate(nCount = strsplit(GeneRatio, split="/") %>% map(~.x[1]) %>% unlist )  %>% 
  mutate(totalCount = strsplit(GeneRatio, split="/") %>% map(~.x[2]) %>% unlist )  %>% 
  dplyr::select(-GeneRatio) %>%
  write.csv(here("results","Saunders","inf_analysis","Kunkle_hc_micro_kegg.csv"),row.names = F,quote=F)

tau_hc.ec_dfbetas %>% write.table(here("results","Saunders","inf_analysis","tau_hc_ec_dfbetas.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ec_go %>% as.data.frame() %>% write.table(here("results","Saunders","inf_analysis","tau_hc_ec_go.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ec_kegg.df %>%  
  ungroup() %>%
  mutate(nCount = strsplit(GeneRatio, split="/") %>% map(~.x[1]) %>% unlist )  %>% 
  mutate(totalCount = strsplit(GeneRatio, split="/") %>% map(~.x[2]) %>% unlist )  %>% 
  dplyr::select(-GeneRatio) %>%
  write.csv(here("results","Saunders","inf_analysis","tau_hc_ec_kegg.csv"),row.names = F,quote=F)

tau_hc.ec2_dfbetas %>% write.table(here("results","Saunders","inf_analysis","tau_hc_ec2_dfbetas.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ec2_go %>% as.data.frame() %>% write.table(here("results","Saunders","inf_analysis","tau_hc_ec2_go.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ec2_kegg.df %>%  
  ungroup() %>%
  mutate(nCount = strsplit(GeneRatio, split="/") %>% map(~.x[1]) %>% unlist )  %>% 
  mutate(totalCount = strsplit(GeneRatio, split="/") %>% map(~.x[2]) %>% unlist )  %>% 
  dplyr::select(-GeneRatio) %>%
  write.csv(here("results","Saunders","inf_analysis","tau_hc_ec2_kegg.csv"),row.names = F,quote=F)

tau_hc.ca1_dfbetas %>% write.table(here("results","Saunders","inf_analysis","tau_hc_ca1_dfbetas.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ca1_go %>% as.data.frame() %>% write.table(here("results","Saunders","inf_analysis","tau_hc_ca1_go.txt"),row.names = F,quote=F, sep="\t")
tau_hc.ca1_kegg.df %>%  
  ungroup() %>%
  mutate(nCount = strsplit(GeneRatio, split="/") %>% map(~.x[1]) %>% unlist )  %>% 
  mutate(totalCount = strsplit(GeneRatio, split="/") %>% map(~.x[2]) %>% unlist )  %>% 
  dplyr::select(-GeneRatio) %>%
  write.csv(here("results","Saunders","inf_analysis","tau_hc_ca1_kegg.csv"),row.names = F,quote=F)

