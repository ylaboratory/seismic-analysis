#null simulation: take input data set
#1: expression data sets
#2: gene set directory. By default 
#3: output file directory
args = commandArgs(trailingOnly = TRUE)
options(warn = -1)

#load packages
suppressMessages(library("here"))
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))

#argument
data_path = args[1]
gs_dir = args[2]
output_file = args[3]

#start processing
suppressMessages(library("seismicGWAS"))

#import data set
suppressMessages(load(data_path))
data("mmu_hsa_mapping")

start = Sys.time()
ts.sce = cal_stat(data_obj = ts.sce, meta_data = as.data.frame(colData(ts.sce)), group = "cell_ontology_class")
ts.sce = cal_sscore(data_obj = ts.sce) 
ts.sce = trans_mmu_to_hsa_stat(ts.sce, gene_mapping_table=mmu_hsa_mapping, from="hsa_ensembl", to="hsa_entrez")
ts.sce = add_glob_stats(ts.sce, stats = c("det_cell_num","ave_exp_ct","max_exp_ct") ) 
end = Sys.time()

#processing time
processing_time = as.numeric(difftime(end, start, units = "secs"))

#import gene set
gs_file_all = list.files(gs_dir, pattern = ".genes.out")

group_time = list()

for (gs_file in gs_file_all){
  gs_file_path = paste0(gs_dir,"/",gs_file)
  zscore = load_zscore(gs_file_path)
  t0 = Sys.time()
  ts.sce = ct_asso(ts.sce, zscore, gene_filter_setting = "det_cell_num>=10& ave_exp_ct > 0.1& max_exp_ct>0.1")
  t1 = Sys.time()
  
  #final time
  group_time = append(group_time, as.numeric(difftime(t1, t0, units = "secs")))
}

save(group_time, processing_time, file = output_file) 