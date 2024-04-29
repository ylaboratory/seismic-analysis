#generate magma expression data set in podaman container 
#parameters:
#1: path of the sce
#2: path of the gene set directory
#3: output file directory
#4: temporary intermediate file path + header (will be removed later)
args = commandArgs(trailingOnly = TRUE)
options(warn = -1)

#argument
data_path = args[1]
gs_dir = args[2]
output_file = args[3]
tmp_file_header = args[4]

#parameters
magma_path = "bin/magma"

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
assay(ts.sce, "cpm") = scuttle::calculateCPM(ts.sce, assay.type = "decontXcounts")

#calculate specific genes 
suppressMessages(library("seismicGWAS"))

#### data set processing for magma
start = Sys.time()
mean_mat <- calc_ct_mean(ts.sce, assay_name = "cpm", ct_label_col = group)
mean_mat_hsa <- translate_gene_ids(t(mean_mat), from = "mmu_symbol")
end = Sys.time()
processing_time = as.numeric(difftime(end, start, units = "secs"))

#write output
print_magma_fuma_tbl(t(mean_mat_hsa), table_type = "MAGMA", main_table_path = paste0(tmp_file_header,".magma.txt"))

#run magma
magma_raw_file_all <- list.files(gs_dir, pattern = ".genes.raw")

group_time <- list()
for(magma_file in magma_raw_file_all){
  magma_raw_path <- paste0(gs_dir,"/",magma_file)
  tm <- system.time(system(paste("/grain/ql29/podman_file/bin/magma","--gene-results", magma_raw_path,"--set-annot",paste0(tmp_file_header,".magma.txt"),
                                "--out",paste0(tmp_file_header, magma_file,".magma.txt"))))
  group_time  <- append(group_time , tm["elapsed"])
}

system(paste("rm",paste0(tmp_file_header, "*")))
save(group_time, processing_time, file = output_file) 

