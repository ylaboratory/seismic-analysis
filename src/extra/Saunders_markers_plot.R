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
if (!require("scran")){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("scran")
  library("scran")
}  #normalize data

if (!require("seismicGWAS")){
  if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

##### load data ######
load(here("data","expr","Saunders","Saunders_processed.rda"))

##### calculate specificity score ####
region_subclass_sscore <- calc_specificity(sce = brain_sce, ct_label_col = "region_subclass", min_avg_exp_ct = 0.01)
region_class_sscore <- calc_specificity(sce = brain_sce, ct_label_col = "region_class", min_avg_exp_ct = 0.01)

#plot data frame
top_sn_genes <- region_class_sscore[,"SN.Neurons"]  %>% sort() %>% tail(n=50)
top_sn_gaba_genes <- region_subclass_sscore[,"SN.GABAergic"]  %>% sort() %>% tail(n=50)
top_sn_gluta_genes <- region_subclass_sscore[,"SN.Glutamatergic"]  %>% sort() %>% tail(n=50)
top_sn_dopa_genes <- region_subclass_sscore[,"SN.Dopaminergic"]  %>% sort() %>% tail(n=50)

top_sn_gene_df <- tibble(gene_name = names(top_sn_genes), specificity = top_sn_genes, 
                         gene_type = ifelse(gene_name %in% names(top_sn_gaba_genes), "GABA specific", "others")) %>%
  mutate(gene_type = ifelse(gene_name %in% names(top_sn_gluta_genes), "VGLUT specific", gene_type)) %>%
  mutate(gene_type = ifelse(gene_name %in% names(top_sn_dopa_genes), "DA specific", gene_type)) %>%
  mutate(gene_name = factor(gene_name, levels = rev(gene_name))) %>%
  mutate(gene_type = factor(gene_type, levels = c("DA specific","GABA specific", "VGLUT specific",  "others")))

top_sn_gaba_df <- tibble(gene_name = names(top_sn_gaba_genes), specificity = top_sn_gaba_genes, gene_type = "GABA specific") %>%
  mutate(gene_name = factor(gene_name, levels = rev(gene_name))) %>%
  mutate(gene_type = factor(gene_type, levels = c("DA specific","GABA specific", "VGLUT specific",  "others")))

top_sn_gluta_df <- tibble(gene_name = names(top_sn_gluta_genes), specificity = top_sn_gluta_genes, gene_type = "VGLUT specific") %>%
  mutate(gene_name = factor(gene_name, levels = rev(gene_name))) %>%
  mutate(gene_type = factor(gene_type, levels = c("DA specific","GABA specific", "VGLUT specific",  "others")))

top_sn_dopa_df <- tibble(gene_name = names(top_sn_dopa_genes), specificity = top_sn_dopa_genes, gene_type = "DA specific") %>%
  mutate(gene_name = factor(gene_name, levels = rev(gene_name))) %>%
  mutate(gene_type = factor(gene_type, levels = c("DA specific","GABA specific", "VGLUT specific",  "others")))

f1 = ggplot(top_sn_gene_df, aes(x = gene_name, y = specificity, fill = gene_type)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() +
  scale_fill_manual(values = set_names(ggsci::pal_npg()(4), c("DA specific","GABA specific", "VGLUT specific",  "others"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 1/ncol(region_class_sscore), linetype = "dashed") +
  scale_y_continuous(limits = c(0,1)) +
  ggtitle("Top 50 genes ranked by specifcity score for SN neurons") +
  theme(plot.title = element_text(hjust = 0.5))

f2 = ggplot(top_sn_dopa_df, aes(x = gene_name, y = specificity, fill = gene_type)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() +
  scale_fill_manual(values = set_names(ggsci::pal_npg()(4), c("DA specific","GABA specific", "VGLUT specific",  "others"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 1/ncol(region_subclass_sscore), linetype = "dashed") +
  scale_y_continuous(limits = c(0,1)) +
  ggtitle("Top 50 genes ranked by specifcity score for SN.Dopaminergic neurons") +
  theme(plot.title = element_text(hjust = 0.5))

f3 = ggplot(top_sn_gluta_df, aes(x = gene_name, y = specificity, fill = gene_type)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() +
  scale_fill_manual(values = set_names(ggsci::pal_npg()(4), c("DA specific","GABA specific", "VGLUT specific",  "others"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 1/ncol(region_subclass_sscore), linetype = "dashed") +
  scale_y_continuous(limits = c(0,1)) + 
  ggtitle("Top 50 genes ranked by specifcity score for SN.Glutamatergic neurons") +
  theme(plot.title = element_text(hjust = 0.5))

f4 = ggplot(top_sn_gaba_df, aes(x = gene_name, y = specificity, fill = gene_type)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() +
  scale_fill_manual(values = set_names(ggsci::pal_npg()(4), c("DA specific","GABA specific", "VGLUT specific",  "others"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 1/ncol(region_subclass_sscore), linetype = "dashed")+
  scale_y_continuous(limits = c(0,1)) +
  ggtitle("Top 50 genes ranked by specifcity score for SN.GABAergic neurons") +
  theme(plot.title = element_text(hjust = 0.5))


cowplot::plot_grid( f1, f2, f4, f3, labels = c("A", "B", "C", "D"), ncol = 1)


