#Analysis of Tabula sapiens data set

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
#load(here("data","expr","Tabula_sapiens","TS_clean.rda"))
load(here("data","expr","Tabula_sapiens","TS_clean.new.rda"))

##### 2. quality control, normalization and label cleaning#######
#num umi has been controlled
ts_obj$mt_counts <- assay(ts_obj, "counts")[grep(pattern = "^MT-", rowData(ts_obj)$gene_symbol), ] %>% colSums()
ts_obj$mt_ratio <- ts_obj$mt_counts/ts_obj$n_counts_UMIs 
ts_obj <- ts_obj[, ts_obj$mt_ratio<=0.1]

###filter genes
ts_obj <- ts_obj[(rowSums( assay(ts_obj, "counts")>0) >=10 &  rowSums( assay(ts_obj, "counts"))>=20),]

#normalization
#clustering
cluster = quickCluster(ts_obj,  assay.type = "counts") 
#calculating normalization factors
size.factor = calculateSumFactors(ts_obj, cluster=cluster, min.mean=0.1,assay.type="counts")
#normalize
ts_obj = logNormCounts(ts_obj, size.factors = size.factor,assay.type = "counts")


##### 3. Calculate specificity score and perform cell type association#######
###specificity score and enrichment
ts_sscore  <- calc_specificity(sce = ts_obj , ct_label_col = "cluster_name", min_avg_exp_ct = 0.01) #for droplet

#map to human genes
ts_sscore_hsa <- translate_gene_ids(ts_sscore, from = "hsa_ensembl")

##enrichment
#magma zscore file
magma_zscore_file <- list.files(here("data","gwas","tm_gwas","zscore"), full.names = T) %>%
  set_names(str_extract(string = ., pattern = "(?<=/)[^/]+$") %>% gsub(pattern = "\\.genes\\.out", replacement = "", x=.))

#enrichment
ts_association <- magma_zscore_file %>%
  map(~get_ct_trait_associations(sscore =ts_sscore_hsa, magma = .x))

##save results
ts_res <- ts_association %>%
  map(~select(.x, cell_type, pvalue)) %>%
  map2(names(.), ~set_colnames(.x, c("cell_type", .y))) %>%
  purrr::reduce(~left_join(.x, .y, by="cell_type"))

write.table(ts_res, here("results","Tabula_sapiens","seismic","ts_res.new.txt"),quote=F, sep="\t", row.names = F)

##save objects for later 
save(ts_obj, file=here("data","expr","Tabula_sapiens","TS_processed.new.rda"))

##save cell ontology mapping
colData(ts_obj) %>% 
  as_tibble() %>%
  distinct(cell_ontology_id, organ_tissue, cluster_name) %>%
  set_colnames(c("cell_ontology_id","tissue","cluster_name")) %>%
  write.table(here("results","Tabula_sapiens","ts_ontology.new.txt"), sep="\t",col.names = T, quote = F)

