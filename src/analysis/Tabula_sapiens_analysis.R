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
load(here("data","expr","Tabula_sapiens","TS_clean.rda"))

##### 2. quality control, normalization and label cleaning#######
#num umi has been controlled
ts_obj$mt_counts = assay(ts_obj, "counts")[grep(pattern = "^MT-", rowData(ts_obj)$gene_symbol), ] %>% colSums()
ts_obj$mt_ratio = ts_obj$mt_counts/ts_obj$n_counts_UMIs 
ts_obj = ts_obj[, ts_obj$mt_ratio<=0.1]

###filter genes
ts_obj = ts_obj[(rowSums( assay(ts_obj, "counts")>0) >=10 &  rowSums( assay(ts_obj, "counts"))>=20),]

#normalization
#clustering
cluster = quickCluster(ts_obj,  assay.type = "counts") 
#calculating normalization factors
size.factor = calculateSumFactors(ts_obj, cluster=cluster, min.mean=0.1,assay.type="counts")
#normalize
ts_obj = logNormCounts(ts_obj, size.factors = size.factor,assay.type = "counts")


##### 3. Calculate specificity score and perform cell type association#######
###specificity score and enrichment without out groups
ts_obj = cal_stat(data_obj = ts_obj , meta_data = as.data.frame(colData(ts_obj)), group = "cluster_name") #cluster_name = organ/tissue + cell ontology

#calculate specificity score
ts_obj  = cal_sscore(data_obj = ts_obj) 

#map to human genes
data("mmu_hsa_mapping")  
ts_obj = trans_mmu_to_hsa_stat(ts_obj , gene_mapping_table=mmu_hsa_mapping, from="hsa_ensembl", to="hsa_entrez")

#add global statistics
ts_obj = add_glob_stats(ts_obj, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 

##enrichment
gwas_zscore = load_zscore(here("data","gwas","tm_gwas","zscore"))
ts_obj = cal_ct_asso(ts_obj, gwas_zscore, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.01& max_exp_ct>0.1")

#save results
ts_res = get_ct_asso(ts_obj, trait_name = "all", asso_model = "linear", merge_output = T)

write.table(ts_res, here("results","Tabula_sapiens","seismic","ts_res.txt"),quote=F, sep="\t", row.names = F)

##save objects for later 
save(ts_obj, file=here("data","expr","Tabula_sapiens","TS_processed.rda"))
