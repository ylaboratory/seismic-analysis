#Analysis of Tabula muris data set
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
library("SingleCellExperiment")
library("EWCE")
library("MAGMA.Celltyping")


### load data
#load the object that's filtered already in the Tabula_muris_analysis scripts
load(here("data","expr","Tabula_muris","TM_processed.rda"))
source(here("src","tools","magma_fuma_file_prep.R"))
source(here("src","tools","sparse_mat_util.R"))

### prepare ewce files
assay(facs_obj_sce_munge, "cpm")  <- scuttle::calculateCPM(facs_obj_sce_munge, assay.type = "counts", size.factors = colSums(assay(facs_obj_sce, "counts")))

##remove rare cell types
sscore <- seismicGWAS::calc_specificity(facs_obj_sce, assay_name = "cpm", ct_label_col = "cluster_name")
facs_obj_sce <- facs_obj_sce[, facs_obj_sce$cluster_name %in% colnames(sscore)]

###option 1: munge sce mat 
#munge matrix
source(here("src","tools","munge_sce_mat.R"))
load(here("data","ref","mapping","mmu_hsa_mapping.rda"))
facs_obj_sce_munge = munge_sce_mat(facs_obj_sce, mapping_df = mmu_hsa_mapping %>% distinct(mmu_symbol, hsa_symbol) %>% drop_na(), assay_name = "cpm")

#ctd files
ctd_human <- generate_celltype_data(exp = assay(facs_obj_sce_munge, "cpm"),
                                   annotLevels = list(cluster_name=facs_obj_sce_munge$cluster_name),
                                   groupName = "cluster_name",
                                   no_cores = 10,
                                   input_species = "human",
                                   return_ctd = T) 
ctd_mouse_human <- generate_celltype_data(exp = assay(facs_obj_sce, "cpm"),
                                         annotLevels = list(cluster_name=facs_obj_sce$cluster_name),
                                         groupName = "cluster_name",
                                         no_cores = 10,
                                         input_species = "mouse",
                                         return_ctd = T) 

ctd_human = ctd_human$ctd
ctd_mouse_human = ctd_mouse_human$ctd

names(ctd_human) = c("cluster_name")
ctd_human_new = map_specificity_to_entrez(ctd_human, annotLevel = "cluster_name",ctd_species = "human")


#check specificity difference
mean_ewce <- ctd_human$ctd[[1]]$mean_exp %>% t()
specificity_ewce <- ctd_human$ctd[[1]]$specificity %>% t()

mean_ours <- calc_ct_mean(facs_obj_sce_munge, assay_name = "cpm", ct_label_col = "cluster_name") %>%
  .[match(rownames(mean_ewce), rownames(.)), ]
specificity_ours <- mean_ours %>% sweep( MARGIN = 2, STATS = Matrix::colSums(.), FUN="/") 

###option 2: drop without uniq mappings

MAGMA_resultsres <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs = here("data","gwas","tm_gwas","tm_gwas_magmact","AF"),
  ctd = ctd_human,
  ctd_species = "human", 
  ctd_name = "cluster_name", 
  run_linear = TRUE, 
  run_top10 = TRUE)

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/grain/ql29/seismic-analysis/bin/magma", sep = ":"))

res = calculate_celltype_associations(ctd = ctd,
                                      magma_dir = here("data","gwas","tm_gwas","tm_gwas_magmact","AF"))
###or firstly munge gene names

