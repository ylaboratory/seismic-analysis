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
mmu_hsa_mapping <- args[5]

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

#load mapping
load(mmu_hsa_mapping)
mmu_hsa_mapping <- mmu_hsa_mapping %>% 
  distinct(hsa_ensembl, hsa_entrez) %>%
  drop_na() %>%
  group_by(hsa_ensembl) %>% 
  filter(n()==1) %>%
  group_by(hsa_entrez) %>%
  filter(n()==1) %>%
  ungroup() %>%
  mutate(hsa_entrez = as.character(hsa_entrez))

#load data set
load(data_path)
#transform to cpm
assay(ts.sce, "cpm") = scuttle::calculateCPM(ts.sce, assay.type = "decontXcounts", size.factors = colSums(assay(ts.sce, "decontXcounts")))

#### data set processing for magma
start = Sys.time()
#calculate mean expression
mean_mat <- calc_ct_mean(ts.sce, assay_name = "decontXcounts", ct_label_col = "cell_ontology_class")
#print("step1",as.numeric(difftime(Sys.time(), start, units = "secs")))

#remove not 1 to 1 gene mapping
mean_mat = mean_mat[, colnames(mean_mat) %in%mmu_hsa_mapping$hsa_ensembl] %>% 
  set_colnames(mmu_hsa_mapping$hsa_entrez[match( colnames(.), mmu_hsa_mapping$hsa_ensembl)])
#print("step2",as.numeric(difftime(Sys.time(), start, units = "secs")))

#filter non-expressed genes
mean_mat = mean_mat[, which(colSums(mean_mat)>0)]
#print("step3",as.numeric(difftime(Sys.time(), start, units = "secs")))

end = Sys.time()
processing_time = as.numeric(difftime(end, start, units = "secs"))

#write output
print_magma_fuma_tbl(mean_mat, table_type = "MAGMA", main_table_path = paste0(tmp_file_header,".magma.txt"))


#run magma
magma_raw_file_all <- list.files(gs_dir, pattern = ".genes.raw")

group_time <- list()
for(magma_file in magma_raw_file_all){
  magma_raw_path <- paste0(gs_dir,"/",magma_file)
  tm <- system.time(system(paste("/grain/ql29/podman_file/bin/magma","--gene-results", magma_raw_path,"--set-annot",paste0(tmp_file_header,".magma.txt"),
                                "--out",paste0(tmp_file_header, magma_file,".magma.txt"))))
  group_time  <- append(group_time , tm["elapsed"])
}

#system(paste("rm",paste0(tmp_file_header, "*")))
save(group_time, processing_time, file = output_file) 

