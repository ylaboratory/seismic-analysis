#null simulation: take input data set
#1: expression data sets
#2: gene set directory. By default 
#3: output file directory
args <- commandArgs(trailingOnly = TRUE)
options(warn = -1)

#load packages
suppressMessages(library("here"))
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))

#argument
data_path <- args[1]
gs_dir <- args[2]
output_file <- args[3]

#start processing
suppressMessages(library("seismicGWAS"))

#import data set
suppressMessages(load(data_path))

start <- Sys.time()
seismic_sscore  <- calc_specificity(sce = ts.sce, ct_label_col = "cell_ontology_class")
seismic_sscore_hsa <- translate_gene_ids(seismic_sscore, from = "hsa_ensembl")
end <- Sys.time()

#processing time
processing_time <- as.numeric(difftime(end, start, units = "secs"))

#import gene set
gs_file_all <- list.files(gs_dir, pattern = ".genes.out")

group_time = list()

for (gs_file in gs_file_all) {
  gs_file_path <- paste0(gs_dir,"/",gs_file)
  zscore <- read.table(gs_file_path, header = T)
  t0 = Sys.time()
  p_value_df <- get_ct_trait_associations(sscore = seismic_sscore_hsa, magma = zscore)
  t1 = Sys.time()
  
  #final time
  group_time <- append(group_time, as.numeric(difftime(t1, t0, units = "secs")))
}

save(group_time, processing_time, file = output_file) 