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
load(here("data","expr","Saunders","Saunders_clean.rda"))

##### 2. quality control, normalization and label cleaning#######
colData(brain_sce)$tot_counts = colSums(assay(brain_sce,"counts")) #total counts
colData(brain_sce)$mito_ratio = colSums(assay(brain_sce,"counts")[which(grepl(x=rownames(brain_sce),pattern = "^mt-")), ])/brain_sce$tot_counts #mitochondria RNA ratio

#filter #genes have been filtered
brain_sce = brain_sce[,brain_sce$tot_counts>=1000 & brain_sce$mito_ratio<=0.1]

#normalization
#cluster
cluster.brain = quickCluster(brain_sce, assay.type = "counts") 
#calculate normalization factors
brain.factor = calculateSumFactors(brain_sce, cluster=cluster.brain, min.mean=0.1)
#normalize
brain_sce = logNormCounts(brain_sce, size.factors = brain.factor )

#create objects for other granularities
brain_sce_subclass = brain_sce
brain_sce_region_subclass = brain_sce
brain_sce_region_class = brain_sce
brain_sce_region_cluster = brain_sce

#add analysis granularity column
brain_sce_region_subclass$region_subclass = paste0(brain_sce_region_subclass$region,".",brain_sce_region_subclass$subclass)
brain_sce_region_class$region_class = paste0(brain_sce_region_class$region,".",brain_sce_region_class$class)
brain_sce_region_cluster$region_cluster = paste0(brain_sce_region_cluster$region,".",brain_sce_region_cluster$cluster_anno)

##### 3. Calculate specificity score and perform cell type association#######
###specificity score for all granularities
brain_sce  = cal_stat(data_obj = brain_sce , meta_data = as.data.frame(colData(brain_sce)), group = "fine_cluster")
brain_sce_subclass = cal_stat(data_obj = brain_sce_subclass , meta_data = as.data.frame(colData(brain_sce_subclass)), group = "subclass")
brain_sce_region_subclass = cal_stat(data_obj = brain_sce_region_subclass , meta_data = as.data.frame(colData(brain_sce_region_subclass)), group = "region_subclass")
brain_sce_region_class = cal_stat(data_obj = brain_sce_region_class , meta_data = as.data.frame(colData(brain_sce_region_class)), group = "region_class")
brain_sce_region_cluster = cal_stat(data_obj = brain_sce_region_cluster , meta_data = as.data.frame(colData(brain_sce_region_cluster)), group = "region_cluster")


#calculagte specificity score
brain_sce  = cal_sscore(data_obj = brain_sce) 
brain_sce_subclass = cal_sscore(data_obj = brain_sce_subclass)
brain_sce_region_subclass = cal_sscore(data_obj = brain_sce_region_subclass)
brain_sce_region_class = cal_sscore(data_obj = brain_sce_region_class)
brain_sce_region_cluster = cal_sscore(data_obj = brain_sce_region_cluster)

#map to human genes
data("mmu_hsa_mapping")  
brain_sce = trans_mmu_to_hsa_stat(brain_sce , gene_mapping_table=mmu_hsa_mapping, from="mmu_symbol", to="hsa_entrez")
brain_sce_subclass = trans_mmu_to_hsa_stat(brain_sce_subclass , gene_mapping_table=mmu_hsa_mapping, from="mmu_symbol", to="hsa_entrez")
brain_sce_region_subclass = trans_mmu_to_hsa_stat(brain_sce_region_subclass , gene_mapping_table=mmu_hsa_mapping, from="mmu_symbol", to="hsa_entrez")
brain_sce_region_class = trans_mmu_to_hsa_stat(brain_sce_region_class , gene_mapping_table=mmu_hsa_mapping, from="mmu_symbol", to="hsa_entrez")
brain_sce_region_cluster = trans_mmu_to_hsa_stat(brain_sce_region_cluster , gene_mapping_table=mmu_hsa_mapping, from="mmu_symbol", to="hsa_entrez")

#add global statistics
brain_sce = add_glob_stats(brain_sce, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 
brain_sce_subclass = add_glob_stats(brain_sce_subclass, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 
brain_sce_region_subclass = add_glob_stats(brain_sce_region_subclass, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 
brain_sce_region_class = add_glob_stats(brain_sce_region_class, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 
brain_sce_region_cluster = add_glob_stats(brain_sce_region_cluster, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 


##enrichment
#pd case
gwas_zscore_pd = load_zscore(here("data","gwas","neuron","zscore","PD.genes.out"))

#enrichment for PD
brain_sce = cal_ct_asso(brain_sce, gwas_zscore_pd, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.01& max_exp_ct>0.1")
brain_sce_subclass = cal_ct_asso(brain_sce_subclass, gwas_zscore_pd, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.01& max_exp_ct>0.1")
brain_sce_region_subclass = cal_ct_asso(brain_sce_region_subclass, gwas_zscore_pd, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.01& max_exp_ct>0.1")
brain_sce_region_class = cal_ct_asso(brain_sce_region_class, gwas_zscore_pd, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.01& max_exp_ct>0.1")
brain_sce_region_cluster = cal_ct_asso(brain_sce_region_cluster, gwas_zscore_pd, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.01& max_exp_ct>0.1")

#Kunkle GWAS and CSF tau
gwas_zscore_AD = load_zscore( here("data","gwas","neuron","zscore",c("tau.genes.out","Kunkleetal_2019.genes.out")))
brain_sce = cal_ct_asso(brain_sce, gwas_zscore_AD , gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.01& max_exp_ct>0.1")

##### 4. export and save the final results ####
##save results
fine_cluster_res = get_ct_asso(brain_sce, trait_name = "all", asso_model = "linear", merge_output = T)
subclass_res = get_ct_asso(brain_sce_subclass, trait_name = "all", asso_model = "linear",merge_output = T)
region_subclass_res = get_ct_asso(brain_sce_region_subclass, trait_name = "all", asso_model = "linear",merge_output = T)
region_classs_res = get_ct_asso(brain_sce_region_class , trait_name = "all", asso_model = "linear",merge_output = T)
region_cluster_res = get_ct_asso(brain_sce_region_cluster , trait_name = "all", asso_model = "linear",merge_output = T)

write.table(fine_cluster_res , here("results","Saunders","fine_cluster","seismic","fine_cluster_res.txt"),quote=F, sep="\t", row.names = F)
write.table(subclass_res, here("results","Saunders","subclass","seismic","subclass_res.txt"),quote=F, sep="\t", row.names = F)
write.table(region_subclass_res , here("results","Saunders","region_subclass","seismic","region_subclass_res.txt"),quote=F, sep="\t", row.names = F)
write.table(region_classs_res , here("results","Saunders","region_class","seismic","region_classs_res.txt"),quote=F, sep="\t", row.names = F)
write.table(region_cluster_res, here("results","Saunders","region_cluster","seismic","region_cluster_res.txt"),quote=F, sep="\t", row.names = F)

#save data for later file preparation
save(brain_sce, brain_sce_subclass, brain_sce_region_subclass, brain_sce_region_class, brain_sce_region_cluster, file = here("data","expr","Saunders","Saunders_processed_multi.rda") )
save(brain_sce, brain_sce_subclass, brain_sce_region_subclass, brain_sce_region_class, brain_sce_region_cluster, file = here("data","expr","Saunders","Saunders_processed.rda") )
