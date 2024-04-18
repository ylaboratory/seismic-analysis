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
facs_obj_sce = facs_obj_sce[,facs_obj_sce$nReads>=2000]
droplet_obj_sce = droplet_obj_sce[,droplet_obj_sce$nCount_RNA>=2000]

#filter cells without cell ontology id
facs_obj_sce = facs_obj_sce[,!is.na(facs_obj_sce$cell_ontology_id)]
droplet_obj_sce = droplet_obj_sce[,!is.na(droplet_obj_sce$cell_ontology_id)]

#filter genes (reduce computational cost)
rowData(facs_obj_sce)$num_cells = rowSums(assay(facs_obj_sce,"counts")>0)
rowData(droplet_obj_sce)$num_cells = rowSums(assay(droplet_obj_sce,"counts")>0)
rowData(facs_obj_sce)$num_counts = rowSums(assay(facs_obj_sce,"counts"))
rowData(droplet_obj_sce)$num_counts = rowSums(assay(droplet_obj_sce,"counts"))

facs_obj_sce = facs_obj_sce[rowData(facs_obj_sce)$num_cells>=5 & rowData(facs_obj_sce)$num_counts>=10,]
droplet_obj_sce = droplet_obj_sce[rowData(droplet_obj_sce)$num_cells>=5 & rowData(droplet_obj_sce)$num_counts>=10,]

#set analysis granularity
#tissue + cell ontology / tissue + free annotation
colData(facs_obj_sce) = colData(facs_obj_sce) %>% 
  as_tibble() %>% #easy to process
  mutate(cluster_name = ifelse(!is.na(free_annotation), paste0(tissue,".",free_annotation), paste0(tissue,".",cell_ontology_class))) %>%
  DataFrame()

colData(droplet_obj_sce) = colData(droplet_obj_sce) %>% 
  as_tibble() %>% #easy to process
  mutate(cluster_name = ifelse(!is.na(free_annotation), paste0(tissue,".",free_annotation), paste0(tissue,".",cell_ontology_class))) %>%
  DataFrame()

#normalization
#clustering
cluster.facs = quickCluster(facs_obj_sce, assay.type = "counts") 
cluster.droplet = quickCluster(droplet_obj_sce, assay.type = "counts")
#calculating normalization factors
facs.factor = calculateSumFactors(facs_obj_sce, cluster=cluster.facs, min.mean=0.1)
droplet.factor = calculateSumFactors(droplet_obj_sce, cluster=cluster.droplet, min.mean=0.1 )
#normalize
facs_obj_sce = logNormCounts(facs_obj_sce, size.factors = facs.factor )
droplet_obj_sce = logNormCounts(droplet_obj_sce, size.factors = droplet.factor )

##### 3. seismic analysis: calculate specificity score and perform cell type association#######
###specificity score and enrichment without out groups
facs_obj_sce  = cal_stat(data_obj = facs_obj_sce , meta_data = as.data.frame(colData(facs_obj_sce)), group = "cluster_name")
droplet_obj_sce  = cal_stat(data_obj = droplet_obj_sce , meta_data = as.data.frame(colData(droplet_obj_sce)), group = "cluster_name")

#calculagte specificity score
facs_obj_sce  = cal_sscore(data_obj = facs_obj_sce) 
droplet_obj_sce  = cal_sscore(data_obj = droplet_obj_sce)

#map to human genes
data("mmu_hsa_mapping")  
facs_obj_sce = trans_mmu_to_hsa_stat(facs_obj_sce , gene_mapping_table=mmu_hsa_mapping, from="mmu_symbol", to="hsa_entrez")
droplet_obj_sce = trans_mmu_to_hsa_stat(droplet_obj_sce , gene_mapping_table=mmu_hsa_mapping, from="mmu_symbol", to="hsa_entrez")

#add global statistics
facs_obj_sce = add_glob_stats(facs_obj_sce, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 
droplet_obj_sce = add_glob_stats(droplet_obj_sce, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 

##enrichment
#load gwas zscore
gwas_zscore = load_zscore(here("data","gwas","tm_gwas","zscore"))

#enrichment
facs_obj_sce = cal_ct_asso(facs_obj_sce, gwas_zscore, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.1& max_exp_ct>0.1", asso_model = "linear")
droplet_obj_sce = cal_ct_asso(droplet_obj_sce, gwas_zscore, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.01& max_exp_ct>0.1",  asso_model = "linear")

##save results
facs_res = get_ct_asso(facs_obj_sce, trait_name = "all", asso_model = "linear", merge_output = T)
droplet_res = get_ct_asso(droplet_obj_sce, trait_name = "all", asso_model = "linear", merge_output = T)

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
