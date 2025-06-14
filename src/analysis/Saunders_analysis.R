# --------------------------------------
# Analysis of Saunders et al dataset
# --------------------------------------

# load packages
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

if (!require("scran")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("scran")
  library("scran")
}

if (!require("seismicGWAS")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

if (!require("anndata")) {
  install.packages("anndata")
  library("anndata")
}

if (!require("clusterProfiler")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("clusterProfiler")
  library("clusterProfiler")
} 

# update path to downloaded files
brain_sce <- readRDS(here("all_data", "processed_expression", "Saunders_processed.rds"))

# the following lines have already been calculated for the processed data
# uncomment and run these when starting with raw data

# # quality control, normalization and label cleaning
# colData(brain_sce)$tot_counts <- colSums(assay(brain_sce, "counts"))
# mito_rows <- grepl(x = rownames(brain_sce), pattern = "^mt-")
# mito_counts <- colSums(assay(brain_sce, "counts")[mito_rows, ])
# colData(brain_sce)$mito_ratio <- mito_counts / brain_sce$tot_count

# # filter genes
# filtered_cells <- brain_sce$tot_counts >= 1000 & brain_sce$mito_ratio <= 0.1
# brain_sce <- brain_sce[, filtered_cells]

# # normalization
# cluster_brain <- quickCluster(brain_sce, assay.type = "counts")
# brain_factor <- calculateSumFactors(brain_sce,
#                                     cluster = cluster_brain,
#                                     min.mean = 0.1)
# brain_sce <- logNormCounts(brain_sce, size.factors = brain_factor)

# # add analysis granularity column
# brain_sce$region_subclass <-
#   paste0(brain_sce$region, ".", brain_sce$subclass)
# brain_sce$region_class <-
#   paste0(brain_sce$region, ".", brain_sce$class)
# brain_sce$region_cluster <-
#   paste0(brain_sce$region, ".", brain_sce$cluster_anno)

# Calculate specificity score and cell type associations
fine_cluster_sscore <- calc_specificity(sce = brain_sce,
                                        ct_label_col = "fine_cluster",
                                        min_avg_exp_ct = 0.01)

subclass_sscore <- calc_specificity(sce = brain_sce,
                                    ct_label_col = "subclass",
                                    min_avg_exp_ct = 0.01)

region_subclass_sscore <- calc_specificity(sce = brain_sce,
                                           ct_label_col = "region_subclass",
                                           min_avg_exp_ct = 0.01)

region_class_sscore <- calc_specificity(sce = brain_sce,
                                        ct_label_col = "region_class",
                                        min_avg_exp_ct = 0.01)

region_cluster_sscore <- calc_specificity(sce = brain_sce,
                                          ct_label_col = "region_cluster",
                                          min_avg_exp_ct = 0.01)

# map to human genes
fine_cluster_sscore_hsa <- translate_gene_ids(fine_cluster_sscore,
                                              from = "mmu_symbol")

subclass_sscore_hsa <- translate_gene_ids(subclass_sscore,
                                          from = "mmu_symbol")

region_subclass_sscore_hsa <- translate_gene_ids(region_subclass_sscore,
                                                 from = "mmu_symbol")

region_class_sscore_hsa <- translate_gene_ids(region_class_sscore,
                                              from = "mmu_symbol")

region_cluster_sscore_hsa <- translate_gene_ids(region_cluster_sscore,
                                                from = "mmu_symbol")

# enrichment calculations
gwas_zscore_pd <- here("all_data", "processed_gwas", "brain_gwas", "zscore", "PD.genes.out")
gwas_zscore_all <- list(
  PD = here("all_data", "processed_gwas", "brain_gwas", "zscore", "PD.genes.out"),
  tau = here("all_data", "processed_gwas", "brain_gwas", "zscore", "tau.genes.out"),
  Kunkleetal_2019 = here("all_data", "processed_gwas", "brain_gwas", "zscore", "Kunkleetal_2019.genes.out")
)


fine_cluster_association <- gwas_zscore_all %>%
  map(~ get_ct_trait_associations(sscore = fine_cluster_sscore_hsa, magma = .x))
subclass_association <- get_ct_trait_associations(sscore = subclass_sscore_hsa, magma = gwas_zscore_pd)

region_subclass_association <- get_ct_trait_associations(sscore = region_subclass_sscore_hsa, magma = gwas_zscore_pd)

region_class_association <- get_ct_trait_associations(sscore = region_class_sscore_hsa, magma = gwas_zscore_pd)

region_cluster_association <- get_ct_trait_associations(sscore = region_cluster_sscore_hsa, magma = gwas_zscore_pd)

# export and save results
fine_cluster_res <- fine_cluster_association %>%
  map(~ select(.x, cell_type, pvalue)) %>%
  map2(names(.), ~ set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~ left_join(.x, .y, by = "cell_type"))

write.table(fine_cluster_res, here("results", "Saunders", "fine_cluster", "seismic", "fine_cluster_res.txt"), quote = F, sep = "\t", row.names = F)

write.table(subclass_association %>% select(cell_type, pvalue) %>% set_colnames(c("cell_type", "PD")),
            here("results", "Saunders", "subclass", "seismic", "subclass_res.txt"),
            quote = F, sep = "\t", row.names = F
)

write.table(region_subclass_association %>% select(cell_type, pvalue) %>% set_colnames(c("cell_type", "PD")),
            here("results", "Saunders", "region_subclass", "seismic", "region_subclass_res.txt"),
            quote = F, sep = "\t", row.names = F
)

write.table(region_class_association %>% select(cell_type, pvalue) %>% set_colnames(c("cell_type", "PD")),
            here("results", "Saunders", "region_class", "seismic", "region_classs_res.txt"),
            quote = F, sep = "\t", row.names = F
)

write.table(region_cluster_association %>% select(cell_type, pvalue) %>% set_colnames(c("cell_type", "PD")),
            here("results", "Saunders", "region_cluster", "seismic", "region_cluster_res.txt"),
            quote = F, sep = "\t", row.names = F
)


# influential gene analysis

# for Kunkle GWAS
Kunkle_res <- fine_cluster_res %>% 
  dplyr::select(cell_type, Kunkleetal_2019) %>% 
  arrange(Kunkleetal_2019) %>% 
  mutate(FDR = p.adjust(Kunkleetal_2019, method="fdr"))

Kunkle_hc.micro_dfbetas <- find_inf_genes("HC.Microglia_Macrophage", fine_cluster_sscore_hsa, here("all_data", "processed_gwas", "brain_gwas", "zscore", "Kunkleetal_2019.genes.out")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))
  
 # for tau GWAS
tau_res <- fine_cluster_res %>% 
  dplyr::select(cell_type, tau) %>% 
  arrange(tau) %>% 
  mutate(FDR = p.adjust(tau, method="fdr"))

tau_hc.ec_dfbetas <- find_inf_genes("HC.Neurons_Entorhinal cortex", fine_cluster_sscore_hsa, here("all_data","processed_gwas","brain_gwas","zscore","tau.genes.out")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))
tau_hc.ec2_dfbetas <- find_inf_genes("HC.Neurons_Medial entorhinal cortex 1", fine_cluster_sscore_hsa, here("all_data","processed_gwas","brain_gwas","zscore","tau.genes.out")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))
tau_hc.ca1_dfbetas <- find_inf_genes("HC.Neurons_CA1 Principal cells", fine_cluster_sscore_hsa, here("all_data","processed_gwas","brain_gwas","zscore","tau.genes.out"))  %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys = gene, keytype = "ENTREZID", column="SYMBOL"))


# for Kunkle GWAS ~ HC.MG
Kunkle_hc.mg_inf_genes <- Kunkle_hc.micro_dfbetas %>%
  filter(zstat>0 & is_influential) %>%
  pull(gene)

Kunkle_hc.mg_go <- enrichGO(gene = Kunkle_hc.mg_inf_genes , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr",pvalueCutoff = 0.2)
Kunkle_hc.mg_kegg <- enrichKEGG(gene = Kunkle_hc.mg_inf_genes, organism = "hsa", pvalueCutoff = 0.2)

Kunkle_hc.mg_kegg.df <- Kunkle_hc.mg_kegg  %>% 
  as.data.frame() %>%
  rowwise() %>%
  mutate(new_geneID = paste0(mapIds(org.Hs.eg.db, keys= unlist(strsplit(geneID, split = "/")), keytype = "ENTREZID", column = "SYMBOL"), collapse = "/"))

# for tau GWAS ~ HC.EC / HC.EC2
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

# save results
# microglia ~ clinical AD
Kunkle_hc.micro_dfbetas %>%
  relocate(gene, genename) %>%
  set_colnames(c("entrez id", "gene symbol", "specificity score", "gene risk (z-score)", "dfbetas", "influential")) %>%
  write.table(here("results", "Saunders", "Kunkle_hc_micro_dfbetas.txt"), row.names = F, quote = F, sep="\t")

Kunkle_hc.mg_go %>% 
  as.data.frame() %>% 
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "GO name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results", "Saunders", "Kunkle_hc_micro_go.txt"), row.names = F, quote = F,sep="\t")

Kunkle_hc.mg_kegg.df %>% 
  ungroup() %>%
  dplyr::select(ID, Description, p.adjust, Count, new_geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "KEGG name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results", "Saunders", "Kunkle_hc_micro_kegg.txt"),row.names = F, quote = F, sep="\t")

# ec ~ tau
tau_hc.ec_dfbetas %>% 
  relocate(gene, genename) %>%
  set_colnames(c("entrez id", "gene symbol", "specificity score", "gene risk (z-score)", "dfbetas", "influential")) %>%
  write.table(here("results", "Saunders", "tau_hc_ec_dfbetas.txt"),row.names = F, quote = F, sep="\t")

tau_hc.ec_go %>% 
  as.data.frame() %>% 
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "GO name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results", "Saunders", "tau_hc_ec_go.txt"),row.names = F, quote = F, sep="\t")

tau_hc.ec_kegg.df %>%  
  ungroup() %>%
  dplyr::select(ID, Description, p.adjust, Count, new_geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "KEGG name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results", "Saunders", "tau_hc_ec_kegg.txt"),row.names = F, quote = F, sep="\t")

# ec2 ~ tau
tau_hc.ec2_dfbetas %>% 
  relocate(gene, genename) %>%
  set_colnames(c("entrez id", "gene symbol", "specificity score", "gene risk (z-score)", "dfbetas", "influential")) %>%
  write.table(here("results", "Saunders", "tau_hc_ec2_dfbetas.txt"),row.names = F, quote = F,  sep="\t")

tau_hc.ec2_go %>% 
  as.data.frame() %>% 
  filter(Count >= 3) %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  set_colnames(c("ID", "GO name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results", "Saunders", "tau_hc_ec2_go.txt"),row.names = F, quote = F, sep="\t")

tau_hc.ec2_kegg.df %>%  
  ungroup() %>%
  dplyr::select(ID, Description, p.adjust,Count, new_geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "KEGG name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results", "Saunders", "tau_hc_ec2_kegg.txt"),row.names = F, quote = F, sep="\t")

# ca1 ~ tau
tau_hc.ca1_dfbetas %>% 
  relocate(gene, genename) %>%
  set_colnames(c("entrez id", "gene symbol", "specificity score", "gene risk (z-score)", "dfbetas", "influential")) %>%
  write.table(here("results", "Saunders", "tau_hc_ca1_dfbetas.txt"), row.names = F, quote = F, sep="\t")

tau_hc.ca1_go %>% 
  as.data.frame() %>% 
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "GO name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results", "Saunders", "tau_hc_ca1_go.txt"), row.names = F, quote = F, sep="\t")

tau_hc.ca1_kegg.df %>%  
  ungroup() %>%
  dplyr::select(ID, Description, p.adjust,Count, new_geneID) %>%
  filter(Count >= 3) %>%
  set_colnames(c("ID", "KEGG name", "FDR", "overlap", "genes enriched")) %>%
  write.table(here("results", "Saunders", "tau_hc_ca1_kegg.txt"), row.names = F, quote = F, sep="\t")

