#Analysis of Tabula muris data set
##### 1. load packages and data#######
###load packages
if (!require("here")) {
  install.packages("here")
  library("here")
}

if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}

if (!require("seismicGWAS")) {
  if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
  }
  devtools::install_github("ylaboratory/seismicGWAS")
  library("seismicGWAS")
}

#color
color_mapping_vec <- c("scDRS" ="#46B8DAFF" , "seismic" = "#D43F3AFF" , "S-MAGMA"="#EEA236FF" , "FUMA"="#5CB85CFF")

###### plot for extended ######
standard_para_df_file<- here("data","expr","causal_sim","sc3_standard","perturbed_expr", "parameter_df.txt")

standard_para_df <- read.table(standard_para_df_file,header=T, sep = " ") 

#load data
standard_res_df <- map2(standard_para_df$output_header, standard_para_df$gene_anno_file, ~ {
  #FUMA results
  fuma_res <- read.table(paste0(.x, ".fuma.gsa.out"), header = T) %>% 
    mutate(FDR = p.adjust(P, method="fdr"))
  fuma_anno <- read.table(paste0(.x, ".fuma.aux.txt"), header = T, sep="\t")
  fuma_ct_variable <- fuma_anno$encoded_name[fuma_anno$cell_type == "target_cell_type"]
  fuma_p <- fuma_res$P[which(fuma_res$VARIABLE != fuma_ct_variable)]
  fuma_pos = length(which(fuma_p <= 0.05))
  fuma_neg = length(fuma_p) - fuma_pos
  
  #S-MAGMA results
  magma_res <-  read.table(paste0(.x,".magma.gsa.out"), header = T) %>% 
    mutate(FDR = p.adjust(P, method="fdr"))
  magma_anno <- read.table(paste0(.x,".magma.aux.txt"), header = T, sep="\t")
  magma_ct_variable <- magma_anno$encoded_name[magma_anno$cell_type == "target_cell_type"]
  magma_p <- magma_res$P[which(magma_res$VARIABLE != magma_ct_variable)]
  magma_pos = length(which(magma_p <= 0.05))
  magma_neg = length(magma_p) - magma_pos
  
  #scdrs results 
  scdrs_res <- read.table(paste0(.x,".scdrs_res.csv"), header = T, sep="\t") %>%
    mutate(P = pnorm(assoc_mcz, lower.tail = F)) %>%
    mutate(FDR = p.adjust(P, method="fdr"))
  scdrs_p <-  scdrs_res$P[which(scdrs_res$group!="target_cell_type")]
  scdrs_pos = length(which(scdrs_p <= 0.05))
  scdrs_neg = length(scdrs_p) - scdrs_pos
  
  #seismic results
  seismic_res <- read.table(paste0(.x,".results.txt"), header = T, sep="\t")
  seismic_p <- seismic_res$pvalue[which(seismic_res$cell_type!="target_cell_type")]
  seismic_pos = length(which(seismic_p <= 0.05))
  seismic_neg = length(seismic_p) - seismic_pos
  
  c(fuma_pos, fuma_neg, magma_pos, magma_neg, scdrs_pos, scdrs_neg, seismic_pos, seismic_neg)
}) %>% 
  list_transpose() %>% 
  set_names(c("fuma_pos", "fuma_neg", "magma_pos", "magma_neg", "scdrs_pos", "scdrs_neg", "seismic_pos", "seismic_neg")) %>%
  as_tibble() %>%
  mutate(GWAS = standard_para_df$gs_name)

summary_res_df <- standard_res_df %>% 
  mutate(fuma_pos_ratio = fuma_pos/(fuma_pos + fuma_neg), 
         magma_pos_ratio = magma_pos/(magma_pos + magma_neg),
         scdrs_pos_ratio = scdrs_pos/(scdrs_pos + scdrs_neg),
         seismic_pos_ratio = seismic_pos/(seismic_pos + seismic_neg)) %>%
  group_by(GWAS) %>%
  mutate(mean_fuma_pos_ratio = mean(fuma_pos_ratio),
         sd_fuma_pos_ratio = sd(fuma_pos_ratio),
         mean_magma_pos_ratio = mean(magma_pos_ratio),
         sd_magma_pos_ratio = sd(magma_pos_ratio),
         mean_scdrs_pos_ratio = mean(scdrs_pos_ratio),
         sd_scdrs_pos_ratio = sd(scdrs_pos_ratio),
         mean_seismic_pos_ratio = mean(seismic_pos_ratio),
         sd_seismic_pos_ratio = sd(seismic_pos_ratio)) %>%
  ungroup() %>%
  select(GWAS,starts_with("mean"), starts_with("sd")) %>%
  distinct()


summary_res_df <- summary_res_df %>%
  mutate(GWAS = recode(GWAS, "gs_1" = "Crohn's disease", "gs_2" = "Depression", "gs_3" = "HDL level", "gs_4" = "Monocyte count", "gs_5" = "Chronic kidney disease", 
                       "gs_6" = "LDL level", "gs_7" = "Type I diabetes", "gs_8" = "Lymphocyte count", "gs_9" = "Cardiovascular diseases", "gs_10" = "Type II diabetes")) %>%
  pivot_longer(cols = -GWAS, names_to = c(".value", "method"), names_pattern = "(mean|sd)_(.+)_pos_ratio") %>%
  mutate(method = recode(method, "fuma" = "FUMA", "magma" = "S-MAGMA", "scdrs" = "scDRS")) %>%
  mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA")))


ggplot(summary_res_df, aes(x= GWAS, y = mean, fill = method,  group = method)) +
  geom_bar(stat = "identity", position=position_dodge(),alpha = 0.8, width = 0.7) +
  #geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)), width = .4, linewidth = .3, position = position_dodge(.9)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,           # Rotate labels 90 degrees
      hjust = 1,            # Right-align the text
      vjust = 0.5           # Center vertically
    )
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.5) +
  ylab("ratio of cell types with P-value < 0.05") +
  ggtitle("Proportion of non-perturbed cell types with p-value < 0.05 in causal simulation") +
  theme(plot.title = element_text(hjust = 0.5))
