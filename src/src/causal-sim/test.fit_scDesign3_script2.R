args <- commandArgs(trailingOnly = TRUE)
options(warn = -1)

output_header <- args[1]
num_cores <- round(as.numeric(args[2]))

suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("scDesign3"))
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))

load(output_header)

sce_para <- extract_para(
  sce = sce,
  marginal_list = sce_marginal,
  n_cores = num_cores,
  family_use = "nb",
  new_covariate = sce_data$newCovariate,
  data = sce_data$dat,
  parallelization = "pbmcmapply")

save(sce, sce_marginal, sce_copula, sce_data, sce_para, file = output_header)

# sc_para_list = list()
# for (i in 1:10){
#     print(i)
#     load(paste0( "data/expr/causal_sim/10_ct_scdesign3_model/expr_ds_",i,".new.rda"))
#     sc_para_list[[paste0("expr_ds_",i)]][["copula"]] = sce_copula$copula_list
#     sc_para_list[[paste0("expr_ds_",i)]][["para"]] = sce_para
# }
# save(sc_para_list, file = "data/expr/causal_sim/10_ct_scdesign3_model/sc_para_list.new.rda")
