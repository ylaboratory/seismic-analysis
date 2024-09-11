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
if (!require("scran")){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("scran")
  library("scran")
}  #normalize data

if (!require("seismicGWAS")){
  if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

###load data
load(here("data","expr","Tabula_muris","facs_clean.rda"))
load(here("data","expr","Tabula_muris","droplet_clean.rda"))

##### 2. quality control, normalization and label cleaning#######
#mtRNA has been filtered out
#filter by RNA counts
facs_obj_sce <- facs_obj_sce[,facs_obj_sce$nReads>=2000]
droplet_obj_sce <- droplet_obj_sce[,droplet_obj_sce$nCount_RNA>=2000]

#filter cells without cell ontology id
facs_obj_sce <- facs_obj_sce[,!is.na(facs_obj_sce$cell_ontology_id)]
droplet_obj_sce <- droplet_obj_sce[,!is.na(droplet_obj_sce$cell_ontology_id)]

#filter genes (reduce computational cost)
rowData(facs_obj_sce)$num_cells <- rowSums(assay(facs_obj_sce,"counts")>0)
rowData(droplet_obj_sce)$num_cells <- rowSums(assay(droplet_obj_sce,"counts")>0)
rowData(facs_obj_sce)$num_counts <- rowSums(assay(facs_obj_sce,"counts"))
rowData(droplet_obj_sce)$num_counts <- rowSums(assay(droplet_obj_sce,"counts"))

facs_obj_sce <- facs_obj_sce[rowData(facs_obj_sce)$num_cells>=5 & rowData(facs_obj_sce)$num_counts>=10,]
droplet_obj_sce <- droplet_obj_sce[rowData(droplet_obj_sce)$num_cells>=5 & rowData(droplet_obj_sce)$num_counts>=10,]

#set analysis granularity
#tissue + cell ontology / tissue + free annotation
colData(facs_obj_sce) <- colData(facs_obj_sce) %>% 
  as_tibble() %>% #easy to process
  mutate(cluster_name = ifelse(!is.na(free_annotation), paste0(tissue,".",free_annotation), paste0(tissue,".",cell_ontology_class))) %>%
  DataFrame()

colData(droplet_obj_sce) <- colData(droplet_obj_sce) %>% 
  as_tibble() %>% #easy to process
  mutate(cluster_name = ifelse(!is.na(free_annotation), paste0(tissue,".",free_annotation), paste0(tissue,".",cell_ontology_class))) %>%
  DataFrame()

#normalization
#clustering
cluster.facs <- quickCluster(facs_obj_sce, assay.type = "counts") 
cluster.droplet <- quickCluster(droplet_obj_sce, assay.type = "counts")
#calculating normalization factors
facs.factor <- calculateSumFactors(facs_obj_sce, cluster=cluster.facs, min.mean=0.1)
droplet.factor <- calculateSumFactors(droplet_obj_sce, cluster=cluster.droplet, min.mean=0.1 )
#normalize
facs_obj_sce <- logNormCounts(facs_obj_sce, size.factors = facs.factor )
droplet_obj_sce <- logNormCounts(droplet_obj_sce, size.factors = droplet.factor )

##### 3. seismic analysis: calculate specificity score and perform cell type association#######
###specificity score 
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
facs_association <- magma_zscore_file %>%
  map(~get_ct_trait_associations(sscore = facs_sscore_hsa, magma = .x))
droplet_association <- magma_zscore_file %>%
  map(~get_ct_trait_associations(sscore = droplet_sscore_hsa, magma = .x))

##save results
facs_res <- facs_association %>%
  map(~select(.x, cell_type, pvalue)) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by="cell_type"))

droplet_res <- droplet_association %>%
  map(~select(.x, cell_type, pvalue)) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by="cell_type"))

write.table(facs_res, here("results","Tabula_muris","FACS","seismic","facs_res.txt"),quote=F, sep="\t", row.names = F)
write.table(droplet_res, here("results","Tabula_muris","droplet","seismic","droplet_res.txt"),quote=F, sep="\t", row.names = F)

##save objects for later 
save(facs_obj_sce, droplet_obj_sce, file=here("data","expr","Tabula_muris","TM_processed.rda"))

##save cell ontology mapping
colData(facs_obj_sce) %>% 
  as_tibble() %>%
  distinct(cell_ontology_id, tissue, cluster_name) %>%
  write.table(here("results","Tabula_muris","FACS","facs_ontology.txt"), sep="\t",col.names = T, quote = F)

colData(droplet_obj_sce) %>% 
  as_tibble() %>%
  distinct(cell_ontology_id, tissue, cluster_name) %>%
  write.table(here("results","Tabula_muris","droplet","droplet_ontology.txt"), sep="\t",col.names = T, quote = F)
