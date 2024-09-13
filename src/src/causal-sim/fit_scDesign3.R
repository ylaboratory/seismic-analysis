#fit a scDesign3 model and save it
#automatically use log Library size and cell type as covariate

#load arguments
#1: sce data set
#2: new coldata
#3: output file name
#4: number of cores to use

args <- commandArgs(trailingOnly = TRUE)
options(warn = -1)

sce_file <- args[1]
new_coldata_file <- args[2]
output_header <- args[3]
num_cores <- round(as.numeric(args[4]))

print(num_cores)

#load package
suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("scDesign3"))
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))

#load file
load(sce_file)
if (!tolower(new_coldata_file) %in% c("na","null","no","false")){
  new_coldata <- read.table(new_coldata_file, header = T, sep="\t") %>% mutate(cell_ontology_class = new_cell_ontology_class)
  colData(sce) <- DataFrame(new_coldata)
}


#run scDesign 3 pipeline
sce$library_size <- colSums(assay(sce, "counts"))
colnames(sce) <- sce$cellid

sce_data <- construct_data(sce = sce ,
                           assay_use = "counts",
                           celltype = "cell_ontology_class",
                           pseudotime = NULL,
                           spatial = NULL,
                           other_covariates = "library_size",
                           corr_by = "1")


sce_marginal <- fit_marginal(data = sce_data,
                             predictor = "gene",
                             mu_formula = "cell_ontology_class + offset(log(library_size))",
                             sigma_formula = "1",
                             family_use = "nb",
                             parallelization = "pbmcmapply",
                             n_cores = num_cores,
                             usebam = TRUE)


set.seed(123)
sce_copula <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = sce_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 15,
  input_data = sce_data$dat,
  parallelization = "pbmcmapply")


sce_para <- extract_para(
   sce = sce,
   marginal_list = sce_marginal,
   n_cores = 10,
   family_use = "nb",
   new_covariate = sce_data$newCovariate,
   data = sce_data$dat,
   parallelization = "pbmcmapply")
 
 save(sce, sce_marginal, sce_copula, sce_data, sce_para, file = output_header)