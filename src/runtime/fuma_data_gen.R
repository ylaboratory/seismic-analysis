#generate magma expression data set in podaman container 
#parameters:
#1: path of the sce
#2: path of the gene set directory
#3: output file directory
#4: temporary intermediate file path + header (will be removed later)
args <- commandArgs(trailingOnly = TRUE)
options(warn = -1)

#argument
data_path <- args[1]
gs_dir <- args[2]
output_file <- args[3]
tmp_file_header <- args[4]

#parameters
magma_path <- "bin/magma"

#load packages
suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("seismicGWAS"))

#load source
source("scripts/magma_fuma_file_prep.R")
source("scripts/sparse_mat_util.R")

#load data set
load(data_path)

#transform to cpm
assay(ts.sce,"logcounts") <- sweep_sparse(assay(ts.sce,"decontXcounts"),margin=2, stats=sparseMatrixStats::colSums2(assay(ts.sce,"decontXcounts"))/1e6,fun = "/") %>%
  transform_sparse(fun = function(x) log(x+1))

#processing
start = Sys.time()
mean_mat <- calc_ct_mean(ts.sce, assay_name = "logcounts", ct_label_col = "cell_ontology_class")
mean_mat_hsa <- translate_gene_ids(t(mean_mat), from = "hsa_ensembl")
end = Sys.time()
processing_time = as.numeric(difftime(end, start, units = "secs"))

#write output
print_magma_fuma_tbl(t(mean_mat_hsa), table_type = "FUMA", main_table_path = paste0(tmp_file_header,".fuma.txt"))

#run magma
magma_raw_file_all <- list.files(gs_dir, pattern = ".genes.raw")
group_time <- list()
for(magma_file in magma_raw_file_all) {
  magma_raw_path <- paste0(gs_dir,"/",magma_file)
  tf <- system.time(system(paste(magma_path,"--gene-results",magma_raw_path,"--gene-covar",paste0(tmp_file_header,".fuma.txt"),
                                "--model condition-hide=Average direction=greater --out",paste0(tmp_file_header, magma_file,".fuma.txt"))))
  group_time  <- append(group_time , tf["elapsed"])
}

#system(paste("rm",paste0(tmp_file_header, "*")))
save(group_time, processing_time, file = output_file) 
