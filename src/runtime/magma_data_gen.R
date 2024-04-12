#generate magma expression data set
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
ts.sce = cal_stat(ts.sce, meta_data = as.data.frame(colData(ts.sce)), group = "cell_ontology_class",assay_name ="cpm",mean_only = T )
ts.sce = trans_mmu_to_hsa_stat(ts.sce, gene_mapping_table = mmu_hsa_mapping, from = "hsa_ensembl", to = "hsa_entrez")
#s_mat = metadata(ts.sce)[["group_info"]][["mean_mat"]]
#s_mat = sweep_sparse(s_mat, margin=2, stats = colSums(s_mat),fun="/")
end = Sys.time()
processing_time = as.numeric(difftime(end, start, units = "secs"))

#write output
#s_tbl = t(s_mat) %>% 
#  as.matrix() %>%
#  as_tibble(rownames = "hsa_entrez") %>%
#  pivot_longer(!hsa_entrez, names_to = "cell_ontology_class", values_to = "specificity") %>%
#  group_by(cell_ontology_class) %>%
#  slice_max(specificity, prop=0.1) %>% 
#  summarize( genes = paste(hsa_entrez, collapse = " ")) %>%
#  mutate(cell_ontology_class = paste0("cluster.",1:n()))
#write.table(s_tbl,file=paste0(tmp_file_header,".magma.txt"), col.names = F, row.names = F, sep=" ",quote=F)

#write output
print_magma_fuma_tbl(ts.sce, table_type = "MAGMA", main_table_path = paste0(tmp_file_header,".magma.txt"))

#run magma
magma_raw_file_all = list.files(gs_dir, pattern = ".genes.raw")

group_time = list()
for(magma_file in magma_raw_file_all){
  magma_raw_path = paste0(gs_dir,"/",magma_file)
  tm = system.time(system(paste("/grain/ql29/podman_file/bin/magma","--gene-results", magma_raw_path,"--set-annot",paste0(tmp_file_header,".magma.txt"),
                                "--out",paste0(tmp_file_header, magma_file,".magma.txt"))))
  group_time  = append(group_time , tm["elapsed"])
}

system(paste("rm",paste0(tmp_file_header, "*")))
save(group_time, processing_time, file = output_file) 

