#script to generate data for causal simulation
if (!require("here")){
  install.packages("here")
  library("here")
}
if (!require("magrittr")){
  install.packages("magrittr")
  library("magrittr")
}
if(!require("tidyverse")){
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("SingleCellExperiment")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SingleCellExperiment")
  library("SingleCellExperiment")
} 
if (!require("anndata")){
  install.packages("anndata")
  library("anndata")
}
library("scDesign3")

