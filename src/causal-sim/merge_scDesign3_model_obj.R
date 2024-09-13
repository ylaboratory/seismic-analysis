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
if (!require("scDesign3", quietly = TRUE)){
  if(!requireNamespace("devtools")){
    install.packages("devtools")
  }
  devtools::install_github("SONGDONGYUAN1994/scDesign3")
  library("scDesign3")
}

#load objects for standard fitted scDesign3 model
sc_para_list = list()


#load objects for fitted scDeisn3 model with 3 cell types

