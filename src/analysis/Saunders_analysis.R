#Analysis of Saunders et al data set

##### 1. load packages and data #######
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

if (!require("scran")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("scran")
  library("scran")
} #normalize data

if (!require("seismicGWAS")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

###load data
load(here("data", "expr", "Saunders", "Saunders_clean.rda"))

##### 2. quality control, normalization and label cleaning #######
#total counts
colData(brain_sce)$tot_counts <- colSums(assay(brain_sce, "counts"))
#mitochondria RNA
mito_rows <- grepl(x = rownames(brain_sce), pattern = "^mt-")
mito_counts <- colSums(assay(brain_sce, "counts")[mito_rows, ])
colData(brain_sce)$mito_ratio <- mito_counts / brain_sce$tot_count

#filter #genes have been filtered
filtered_cells <- brain_sce$tot_counts >= 1000 & brain_sce$mito_ratio <= 0.1
brain_sce <- brain_sce[, filtered_cells]

#normalization
#cluster
cluster_brain <- quickCluster(brain_sce, assay.type = "counts")
#calculate normalization factors
brain_factor <- calculateSumFactors(brain_sce,
                                    cluster = cluster_brain,
                                    min.mean = 0.1)
#normalize
brain_sce <- logNormCounts(brain_sce, size.factors = brain_factor)

#add analysis granularity column
brain_sce$region_subclass <-
  paste0(brain_sce$region, ".", brain_sce$subclass)
brain_sce$region_class <-
  paste0(brain_sce$region, ".", brain_sce$class)
brain_sce$region_cluster <-
  paste0(brain_sce$region, ".", brain_sce$cluster_anno)

##### 3. Calculate specificity score and perform cell type association #######
###specificity score for all granularities
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

#map to human genes
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

##enrichment
#pd case
gwas_zscore_pd <- here("data", "gwas", "brain_gwas", "zscore", "PD.genes.out")
gwas_zscore_all <- list(
  PD = here("data", "gwas", "brain_gwas", "zscore", "PD.genes.out"),
  tau = here("data", "gwas", "brain_gwas", "zscore", "tau.genes.out"),
  Kunkleetal_2019 = here("data", "gwas", "brain_gwas", "zscore", "Kunkleetal_2019.genes.out")
)

#enrichment for PD
fine_cluster_association <- gwas_zscore_all %>%
  map(~ get_ct_trait_associations(sscore = fine_cluster_sscore_hsa, magma = .x))
subclass_association <- get_ct_trait_associations(sscore = subclass_sscore_hsa, magma = gwas_zscore_pd)

region_subclass_association <- get_ct_trait_associations(sscore = region_subclass_sscore_hsa, magma = gwas_zscore_pd)

region_class_association <- get_ct_trait_associations(sscore = region_class_sscore_hsa, magma = gwas_zscore_pd)

region_cluster_association <- get_ct_trait_associations(sscore = region_cluster_sscore_hsa, magma = gwas_zscore_pd)

##### 4. export and save the final results ####
##save results
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

#save data for later file preparation
save(brain_sce, file = here("data", "expr", "Saunders", "Saunders_processed.rda"))
