#generate magma expression data set
#1: path of the sce
#2: path of the output 
args = commandArgs(trailingOnly = TRUE)
options(warn = -1)

#argument
data_path = args[1]
gs_dir = args[2]
output_file = args[3]
tmp_file_header = args[4]


#load packages
suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("seismicGWAS"))

#load data set
load(data_path)

#define function
sparse_matrix_transform = function(x, fun = function(x) log2(x)) {
  if (!inherits(x, c("dgTMatrix", "dgCMatrix", "dgRMatrix"))) {
    stop("Only several sparse matrix types are supported: dgTMatrix, dgCMatrix, dgRMatrix")
  }
  # no need to use "match" 
  x@x = fun(x@x)
  return(x)
}

#transform to cpm
#cpm_value = scuttle::calculateCPM(ts.sce, assay.type = "decontXcounts")
#assay(ts.sce,"logcounts") = as(log2(cpm_value+1),"CsparseMatrix")
assay(ts.sce,"logcounts") = sweep_sparse(assay(ts.sce,"decontXcounts"),margin=2, stats=sparseMatrixStats::colSums2(assay(ts.sce,"decontXcounts"))/1e6,fun = "/") %>%
  sparse_matrix_transform(fun = function(x) log(x+1))
#log TPM+1

#processing
start = Sys.time()
ts.sce = cal_stat(ts.sce, meta_data = as.data.frame(colData(ts.sce)), group = "cell_ontology_class" , assay_name = "logcounts")
ts.sce = trans_mmu_to_hsa_stat(ts.sce, gene_mapping_table = mmu_hsa_mapping, from = "hsa_ensembl", to = "hsa_entrez")
mean_mat = metadata(ts.sce)[["group_info"]][["mean_mat"]]
mean_mat = rbind(mean_mat, colMeans(mean_mat))
end = Sys.time()
rownames(mean_mat)[nrow(mean_mat)] = "Average"
processing_time = as.numeric(difftime(end, start, units = "secs"))

#write output
mean_tbl = t(mean_mat) %>% 
  as.matrix() %>%
  as_tibble(rownames = "GENE") %>%
  set_colnames(c("GENE", paste0("cluster.",1:(ncol(.)-2) ), "Average" ))

write.table(mean_tbl,file=paste0(tmp_file_header,".fuma.txt"), col.names =  T, row.names = F, sep=" ",quote=F)

#run magma
magma_raw_file_all = list.files(gs_dir, pattern = ".genes.raw")
group_time = list()
for(magma_file in magma_raw_file_all){
  magma_raw_path = paste0(gs_dir,"/",magma_file)
  tf = system.time(system(paste("/grain/ql29/podman_file/bin/magma","--gene-results",magma_raw_path,"--gene-covar",paste0(tmp_file_header,".fuma.txt"),
                                "--model condition-hide=Average direction=greater --out",paste0(tmp_file_header, magma_file,".fuma.txt"))))
  group_time  = append(group_time , tf["elapsed"])
}

#system(paste("rm",paste0(tmp_file_header, "*")))
save(group_time, processing_time, file = output_file) 
