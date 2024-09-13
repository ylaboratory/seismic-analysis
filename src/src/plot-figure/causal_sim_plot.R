#load results
##### load packages and results ####
if (!require("here")){
  install.packages("here")
  library("here")
}
if (!require("tidyverse")){
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("magrittr")){
  install.packages("magrittr")
  library("magrittr")
}
if (!require("here")){
  install.packages("here")
  library("here")
}

#color
color_mapping_vec <- c("scDRS" ="#46B8DAFF" , "seismic" = "#D43F3AFF" , "S-MAGMA"="#EEA236FF" , "FUMA"="#5CB85CFF")

###### plot for extended ######
extended_para_df_file_1000 <- here("data","expr","causal_sim","sc3_prob_extended","tot_1000","perturbed_expr", "parameter_df.txt")
extended_para_df_file_2000 <- here("data","expr","causal_sim","sc3_prob_extended","tot_2000","perturbed_expr", "parameter_df.txt")

extended_para_df_1000 <- read.table(extended_para_df_file_1000,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
extended_para_df_2000 <- read.table(extended_para_df_file_2000,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 2000)

extended_para_df <- rbind( extended_para_df_1000, extended_para_df_2000)

pr_curve_data <- list()
extended_res_df <- map2(extended_para_df$output_header, extended_para_df$gene_anno_file, ~{
  
  fuma_res <- read.table(paste0(.x, ".fuma.gsa.out"), header = T) %>% mutate(FDR = p.adjust(P, method="fdr"))
  fuma_anno <- read.table(paste0(.x, ".fuma.aux.txt"), header = T, sep="\t")
  fuma_ct_variable <- fuma_anno$encoded_name[fuma_anno$cell_type == "target_cell_type"]
  fuma_p <- fuma_res$P[which(fuma_res$VARIABLE == fuma_ct_variable)]
  fuma_fdr <- fuma_res$FDR[which(fuma_res$VARIABLE == fuma_ct_variable)]
  fuma_ct_rank <- fuma_res %>% arrange(P) %>% mutate(rank = 1:n()) %>% filter(VARIABLE == fuma_ct_variable) %>% pull(rank)
  
  magma_res <-  read.table(paste0(.x,".magma.gsa.out"), header = T) %>% mutate(FDR = p.adjust(P, method="fdr"))
  magma_anno <- read.table(paste0(.x,".magma.aux.txt"), header = T, sep="\t")
  magma_ct_variable <- magma_anno$encoded_name[magma_anno$cell_type == "target_cell_type"]
  magma_p <- magma_res$P[which(magma_res$VARIABLE == magma_ct_variable)]
  magma_fdr <- magma_res$FDR[which(magma_res$VARIABLE == magma_ct_variable)]
  magma_ct_rank <- magma_res %>% arrange(P) %>% mutate(rank = 1:n()) %>% filter(VARIABLE == magma_ct_variable) %>% pull(rank)
  
  seismic_res <- read.table(paste0(.x,".results.txt"), header = T, sep="\t")
  seismic_p <- seismic_res$pvalue[which(seismic_res$cell_type=="target_cell_type")]
  seismic_fdr <- seismic_res$FDR[which(seismic_res$cell_type=="target_cell_type")]
  seismic_ct_rank <- seismic_res %>% arrange(pvalue) %>% mutate(rank = 1:n()) %>% filter(cell_type == "target_cell_type") %>% pull(rank)
  
  gene_anno <- read.table(.y, header = T)
  inf_genes <- read.table(paste0(.x,".inf_genes.txt"), header = T, sep="\t") %>%
    left_join(gene_anno, by=c("gene" = "GENE")) %>% 
    drop_na(is_causal, dfbetas)
  pos_dfbetas <- inf_genes$dfbetas[inf_genes$is_causal]
  neg_dfbetas <- inf_genes$dfbetas[!inf_genes$is_causal]
  seismic_roc <- PRROC::roc.curve(pos_dfbetas, neg_dfbetas)
  seismic_pr <- PRROC::pr.curve(pos_dfbetas, neg_dfbetas, curve = T)
  pr_curve_data[[.x]] <<- seismic_pr$curve
  seismic_auroc <- seismic_roc$auc
  seismic_auprc <- seismic_pr$auc.integral
  

  inf_gene_num <- inf_genes %>% filter(is_influential) %>% nrow()
  causal_inf_gene_num <- inf_genes %>% filter(is_causal, is_influential) %>% nrow()
  tot_gene_num <- inf_genes %>% nrow()
  non_causal_inf_gene_num <- inf_genes %>% filter(!is_causal, !is_influential) %>% nrow()
  #seismic_causal_ratio <- (nrow(filter(inf_genes, is_causal))/ nrow(inf_genes))
  
  # gene_anno <- read.table(.y,header = T)
  # inf_genes <- read.table(paste0(.x,".inf_genes.txt"), header = T, sep="\t") %>%
  #   left_join(gene_anno, by=c("gene" = "GENE")) %>% 
  #   filter(is_influential)
  # seismic_causal_ratio <- (nrow(filter(inf_genes, is_causal))/ nrow(inf_genes))
  
  if (file.exists(paste0(.x,".scdrs_res.csv"))){
    scdrs_res <- read.table(paste0(.x,".scdrs_res.csv"), header = T, sep="\t") %>%
      mutate(P = pnorm(assoc_mcz, lower.tail = F)) %>%
      mutate(FDR = p.adjust(P, method="fdr"))
    scdrs_p <-  scdrs_res$P[which(scdrs_res$group=="target_cell_type")]
    scdrs_fdr <-  scdrs_res$FDR[which(scdrs_res$group=="target_cell_type")]
    scdrs_ct_rank <- scdrs_res %>% arrange(P) %>% mutate(rank = 1:n()) %>% filter(group == "target_cell_type") %>% pull(rank)
  }else{
    scdrs_p <- NA
    scdrs_fdr <- NA
    scdrs_ct_rank <- NA
  }
  c(seismic_fdr, seismic_p, seismic_ct_rank,inf_gene_num, causal_inf_gene_num, tot_gene_num, non_causal_inf_gene_num, 
    seismic_auroc, seismic_auprc,
    scdrs_fdr, scdrs_p,scdrs_ct_rank, fuma_p, fuma_fdr,fuma_ct_rank, magma_p, magma_fdr,magma_ct_rank)
})%>% list_transpose() %>%
  set_names(c("seismic_fdr","seismic_p", "seismic_ct_rank", "seismic_n_inf_gene","seismic_n_causal_inf_gene", "tot_gene_num","seismic_n_non_causal_inf_gene", 
              "seismic_auroc","seismic_auprc","zstat_auprc",
              "scdrs_fdr", "scdrs_p","scdrs_ct_rank", "fuma_p", "fuma_fdr","fuma_ct_rank", "magma_p", "magma_fdr","magma_ct_rank")) %>%
  as_tibble() %>%
  mutate(effect_size = strsplit(extended_para_df$es_name, split="_") %>% map(~.x[2]) %>% unlist) %>%
  mutate(causal_gene_ratio = strsplit(extended_para_df$gr_name, split="_") %>% map(~.x[2]) %>% unlist) %>%
  mutate(data_set = extended_para_df$expr_name, gene_set = extended_para_df$gs_name, gene_num = extended_para_df$gene_num) %>%
  mutate(n_causal_genes = tot_gene_num - seismic_n_inf_gene - seismic_n_non_causal_inf_gene)


extended_power_df <- extended_res_df  %>% 
  group_by(effect_size, causal_gene_ratio, gene_num, gene_set) %>%
  mutate(seismic_ratio = length(which(seismic_fdr<=0.05))/length(which(!is.na(seismic_fdr))), 
         fuma_ratio = length(which(fuma_fdr<=0.05))/length(which(!is.na(fuma_fdr))),
         magma_ratio = length(which(magma_fdr<=0.05))/length(which(!is.na(magma_fdr))),
         scdrs_ratio = length(which(scdrs_fdr<=0.05))/length(which(!is.na(scdrs_fdr)))) %>%
  distinct(gene_num, causal_gene_ratio, effect_size,gene_set, seismic_ratio, fuma_ratio, magma_ratio, scdrs_ratio) %>%
  set_colnames(c("tot_gene_num","causal_gene_ratio","effect_size","gene_set",  "seismic","FUMA", "S-MAGMA","scDRS")) %>%
  pivot_longer(cols = !effect_size & ! causal_gene_ratio & !tot_gene_num &!gene_set, names_to = "method", values_to = "power") %>%
  group_by(effect_size, causal_gene_ratio, tot_gene_num, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  select(-gene_set) %>%
  distinct() %>%
  mutate(causal_gene_num = paste0(1000*as.numeric(causal_gene_ratio)," causal genes"))

# ggplot(extended_power_df, aes(x = effect_size, y = mean_power,fill=method)) +
#   geom_bar(stat = "identity", position=position_dodge()) +
#   facet_grid(rows = vars(causal_gene_ratio), cols = vars(tot_gene_num)) +
#   scale_fill_manual(values = color_mapping_vec) +
#   theme_classic() +
#   ggtitle("Random disease gene sets sampled with Z-score weight") +
#   theme(plot.title = element_text(hjust = 0.5))

#power plot
#figure main 
ggplot(extended_power_df %>% filter(causal_gene_ratio == 0.5, tot_gene_num == 1000) %>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,linewidth=.3, position=position_dodge(.9)) + 
  #facet_grid(rows = vars(causal_gene_ratio), cols = vars(tot_gene_num)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  #ggtitle("Random disease gene sets sampled with Z-score weight") +
  theme(plot.title = element_text(hjust = 0.5))

#figure supplementary
ggplot(extended_power_df %>% filter(causal_gene_ratio %in% c(0.2,0.4,0.5,0.6,0.8), tot_gene_num == 1000) %>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(causal_gene_num)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  #ggtitle("Random disease gene sets sampled with Z-score weight") +
  theme(plot.title = element_text(hjust = 0.5))

#plot for the ranking
#top 1
extended_rank_1_df <- extended_res_df  %>% 
  group_by(effect_size, causal_gene_ratio, gene_num, gene_set) %>%
  mutate(seismic_ratio = length(which(seismic_ct_rank<=1))/length(which(!is.na(seismic_ct_rank))), 
         fuma_ratio = length(which(fuma_ct_rank<=1))/length(which(!is.na(fuma_ct_rank))),
         magma_ratio = length(which(magma_ct_rank<=1))/length(which(!is.na(magma_ct_rank))),
         scdrs_ratio = length(which(scdrs_ct_rank<=1))/length(which(!is.na(scdrs_ct_rank)))) %>%
  distinct(gene_num, causal_gene_ratio, effect_size,gene_set, seismic_ratio, fuma_ratio, magma_ratio, scdrs_ratio) %>%
  set_colnames(c("tot_gene_num","causal_gene_ratio","effect_size","gene_set",  "seismic","FUMA", "S-MAGMA","scDRS")) %>%
  pivot_longer(cols = !effect_size & ! causal_gene_ratio & !tot_gene_num &!gene_set, names_to = "method", values_to = "power") %>%
  group_by(effect_size, causal_gene_ratio, tot_gene_num, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  select(-gene_set) %>%
  distinct()  %>%
  mutate(causal_gene_num = paste0(1000*as.numeric(causal_gene_ratio), " causal genes"))

#top 1 supplementary plot
ggplot(extended_rank_1_df %>% filter(tot_gene_num==1000, causal_gene_ratio %in% c(0.2, 0.4,0.5,0.6,0.8)) %>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(causal_gene_num)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
 #ggtitle("Power of ranking the target cell type as top 1") +
  theme(plot.title = element_text(hjust = 0.5))


#top 5
extended_rank_5_df <- extended_res_df  %>% 
  group_by(effect_size, causal_gene_ratio, gene_num, gene_set) %>%
  mutate(seismic_ratio = length(which(seismic_ct_rank<=5))/length(which(!is.na(seismic_ct_rank))), 
         fuma_ratio = length(which(fuma_ct_rank<=5))/length(which(!is.na(fuma_ct_rank))),
         magma_ratio = length(which(magma_ct_rank<=5))/length(which(!is.na(magma_ct_rank))),
         scdrs_ratio = length(which(scdrs_ct_rank<=5))/length(which(!is.na(scdrs_ct_rank)))) %>%
  distinct(gene_num, causal_gene_ratio, effect_size,gene_set, seismic_ratio, fuma_ratio, magma_ratio, scdrs_ratio) %>%
  set_colnames(c("tot_gene_num","causal_gene_ratio","effect_size","gene_set",  "seismic","FUMA", "S-MAGMA","scDRS")) %>%
  pivot_longer(cols = !effect_size & ! causal_gene_ratio & !tot_gene_num &!gene_set, names_to = "method", values_to = "power") %>%
  group_by(effect_size, causal_gene_ratio, tot_gene_num, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  select(-gene_set) %>%
  distinct()  %>%
  mutate(causal_gene_num = paste0(1000*as.numeric(causal_gene_ratio), " causal genes"))

#top 1 supplementary plot
ggplot(extended_rank_5_df %>% filter(tot_gene_num==1000, causal_gene_ratio %in% c(0.2, 0.4,0.5,0.6,0.8)) %>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(causal_gene_num)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  #ggtitle("Power of ranking the target cell type as top 1") +
  theme(plot.title = element_text(hjust = 0.5))


#plot for the auroc/auprc
# causal_auprc_df <- extended_res_df %>%
#   #filter(seismic_fdr <= 0.05) %>%
#   mutate(n_causal_genes = tot_gene_num - seismic_n_inf_gene - seismic_n_non_causal_inf_gene, random_auprc = n_causal_genes/tot_gene_num, 
#          seismic_auprc = seismic_auprc/random_auprc, zstat_auprc = zstat_auprc/random_auprc) %>%
#   group_by(gene_num, effect_size, causal_gene_ratio) %>% 
#   mutate(mean_auprc = mean(seismic_auprc), sd_auprc= sd(seismic_auprc), mean_z_auprc = mean(zstat_auprc), sd_z_auprc = sd(zstat_auprc)) %>%
#   distinct(gene_num, effect_size, causal_gene_ratio, mean_auprc, sd_auprc, mean_z_auprc, sd_z_auprc) %>%
#   ungroup() 
# 
# causal_auprc_df_long <- causal_auprc_df %>%
#   pivot_longer(cols = starts_with("mean_") | starts_with("sd_"), names_to = c(".value", "auprc_fold_change"), names_pattern = "(mean|sd)_(.*)")
# 
# ggplot(causal_auprc_df_long, aes(x = effect_size, y = log2(mean), group = auprc_fold_change , fill = auprc_fold_change )) +
#   geom_bar(stat="identity", position = position_dodge(),alpha = 0.3) +
#   geom_errorbar(aes(ymin = log2(pmax(0, mean - sd)), ymax = log2(mean + sd)), width=.2, position = position_dodge(0.9)) + 
#   facet_grid(rows = vars(causal_gene_ratio), cols = vars(gene_num)) +
#   ylab("log2 auroc fold change (compared with random baseline)")+
#   theme_minimal() +
#   ggtitle("AUROC for assessing causal genes via DFBETAS") +
#   theme(plot.title = element_text(hjust = 0.5))

#PRC plot: get the condition
pr_curve_data <- pr_curve_data %>% 
  map(~set_colnames(.x, c("recall", "precision","score"))) %>% 
  map(~as_tibble(.x))

plot_pr_curve_data <- pr_curve_data %>% 
  map2(extended_res_df$n_causal_genes / extended_res_df$tot_gene_num, ~mutate(.x, precision = precision/.y)) %>% 
  .[extended_res_df$gene_num == 1000 & extended_res_df$causal_gene_ratio == 0.5] %>% 
  split(extended_res_df$effect_size[extended_res_df$gene_num == 1000 & extended_res_df$causal_gene_ratio == 0.5]) %>%
  map(~map(.x, ~approx(x = .x[["recall"]], y = .x[["precision"]], xout = seq(0.01,1,0.01)))) %>%
  map(~map(.x, ~as_tibble(.x))) %>%
  map(~map(.x, ~set_colnames(.x, c("recall", "precision_fold")))) %>% 
  map(~reduce(.x, ~rbind(.x, .y))) %>%
  map(~group_by(.x, recall)) %>%
  map(~mutate(.x, mean_precision_fold = mean(precision_fold), sd_precision_fold = sd(precision_fold))) %>%
  map(~ungroup(.x)) %>%
  map(~distinct(.x, recall, mean_precision_fold, sd_precision_fold)) %>%
  map2(names(.), ~mutate(.x, effect_size = as.numeric(.y))) %>%
  reduce(~rbind(.x, .y ))

ggplot(plot_pr_curve_data %>% mutate(effect_size = as.character(effect_size)), aes(x = recall, y = log2(mean_precision_fold), color = effect_size, group = effect_size)) +
  geom_line() +
  #geom_ribbon(aes(ymin = log2(pmax(0.001, mean_precision_fold - sd_precision_fold)), 
  #                ymax = log2(mean_precision_fold + sd_precision_fold), 
  #                fill = effect_size), 
  #            alpha = 0.3, color=NA) +
  ylab("precision over random") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")) +
  theme_classic() +
  scale_color_manual(values = rev(viridis::viridis(6, alpha=0.8))) +
  guides(color = guide_legend(title = "log2(expression fold change)"))

#plot for causal ratio
# causal_ratio_df <- extended_res_df %>%
#   mutate(seismic_causal_ratio = seismic_n_causal_inf_gene/seismic_n_inf_gene) %>%
#   group_by(gene_num, effect_size, causal_gene_ratio) %>%
#   #filter(seismic_fdr <= 0.05) %>%
#   mutate(mean_ratio = mean(seismic_causal_ratio), sd_ratio= sd(seismic_causal_ratio)) %>%
#   distinct(gene_num, effect_size, causal_gene_ratio, mean_ratio, sd_ratio) %>%
#   ungroup() 
# 
# ggplot(causal_ratio_df, aes(x = effect_size, y = mean_ratio)) +
#   geom_bar(stat="identity", position = position_dodge(),alpha = 0.3) +
#   geom_errorbar(aes(ymin = pmax(0, mean_ratio - sd_ratio), ymax = pmin(1, mean_ratio + sd_ratio)), width=.1, linewidth = 0.2) + 
#   facet_grid(rows = vars(causal_gene_ratio), cols = vars(gene_num)) +
#   ylab("ratio")+
#   theme_minimal() +
#   ggtitle("Ratio of influential genes that are causal genes") +
#   theme(plot.title = element_text(hjust = 0.5))


 #plot for enrichment fold
enrich_fold_df <- extended_res_df %>%
  mutate(n_causal_genes = tot_gene_num - seismic_n_inf_gene - seismic_n_non_causal_inf_gene, 
         enrichment_ratio = (seismic_n_causal_inf_gene/seismic_n_inf_gene)/(n_causal_genes/tot_gene_num),
         z_enrichment_ratio = (top_z_cuasal_gene_num/seismic_n_inf_gene) / (n_causal_genes/tot_gene_num),
         max_enrichment_ratio = (pmin(n_causal_genes, seismic_n_causal_inf_gene) / seismic_n_causal_inf_gene)/(n_causal_genes/tot_gene_num)) %>%
  select(gene_num, effect_size, causal_gene_ratio, enrichment_ratio,  z_enrichment_ratio, max_enrichment_ratio) %>%
  group_by(gene_num, effect_size, causal_gene_ratio) %>%
  mutate(mean_enrichment_ratio = mean(enrichment_ratio), sd_enrichment_ratio = sd(enrichment_ratio)) %>%
  mutate(mean_z_enrichment_ratio = mean(z_enrichment_ratio), sd_z_enrichment_ratio = sd(z_enrichment_ratio)) %>%
  mutate(mean_max_enrichment_ratio = mean(max_enrichment_ratio), sd_max_enrichment_ratio = sd(max_enrichment_ratio)) %>%
  distinct(gene_num, effect_size, causal_gene_ratio, mean_enrichment_ratio, sd_enrichment_ratio, 
           mean_z_enrichment_ratio, sd_z_enrichment_ratio, mean_max_enrichment_ratio, sd_max_enrichment_ratio) %>%
  ungroup() 

enrich_fold_df_long <- enrich_fold_df %>% 
  pivot_longer(cols = starts_with("mean_") | starts_with("sd_"), names_to = c(".value", "enrichment_type"), names_pattern = "(mean|sd)_(.*)")

ggplot(enrich_fold_df, aes(x = effect_size, y = log2(mean_enrichment_ratio))) +
  geom_bar(stat="identity", position = position_dodge(),alpha = 0.3) +
  geom_errorbar(aes(ymin =  log2(mean_enrichment_ratio - sd_enrichment_ratio), ymax = log2(mean_enrichment_ratio + sd_enrichment_ratio)), width=.2) + 
  facet_grid(rows = vars(causal_gene_ratio), cols = vars(gene_num), scales = "free_y") +
  ylab("Enrichment ratio (log2 Scale)")+
  theme_minimal() +
  ggtitle("Enrichment ratio of causal genes for influential genes compared with others") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(enrich_fold_df_long, aes(x = effect_size, y = log2(mean), group = enrichment_type, fill = enrichment_type)) +
  geom_bar(stat="identity", position = position_dodge(), alpha = 0.6) +
  geom_errorbar(aes(ymin =  log2(mean - sd), ymax = log2(mean + sd)), width=.2, position = position_dodge(0.9),linewidth = 0.3) + 
  facet_grid(rows = vars(causal_gene_ratio), cols = vars(gene_num), scales = "free_y") +
  ylab("Enrichment ratio (log2 Scale)")+
  theme_minimal() +
  ggtitle("Enrichment ratio of causal genes for influential genes compared with others") +
  theme(plot.title = element_text(hjust = 0.5))

# ###### plot for multiple cell types ######
# extended_para_df_file_random <- here("data","expr","causal_sim","sc3_multi_ct","random","perturbed_expr", "parameter_df.txt")
# extended_para_df_file_no_overlap <- here("data","expr","causal_sim","sc3_multi_ct","no_overlap","perturbed_expr", "parameter_df.txt")
# extended_para_df_file_causal_overlap <- here("data","expr","causal_sim","sc3_multi_ct","causal_overlap","perturbed_expr", "parameter_df.txt")
# extended_para_df_file_ct_overlap <- here("data","expr","causal_sim","sc3_multi_ct","ct_overlap","perturbed_expr", "parameter_df.txt")
# 
# extended_para_df_random <- read.table(extended_para_df_file_random,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
# extended_para_df_no_overlap <- read.table(extended_para_df_file_no_overlap,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
# extended_para_df_causal_overlap <- read.table(extended_para_df_file_causal_overlap,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
# extended_para_df_ct_overlap <- read.table(extended_para_df_file_ct_overlap,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
# 
# #extended_para_df <- rbind(extended_para_df_500, extended_para_df_1000, extended_para_df_2000, extended_para_df_5000)
# extended_para_df_multi <- rbind(extended_para_df_random, extended_para_df_no_overlap, extended_para_df_causal_overlap, extended_para_df_ct_overlap)
# 
# extended_res_df_multi <- map2(extended_para_df_multi$output_header, 1:nrow(extended_para_df_multi), ~{
#   
#   if(.y %% 100 == 0){
#     print(.y)
#   }
#   
#   target_cell_type = c("target_cell_type_1", "target_cell_type_2", "target_cell_type_3")
#   
#   fuma_res <- read.table(paste0(.x, ".fuma.gsa.out"), header = T) %>% mutate(FDR = p.adjust(P, method="fdr"))
#   fuma_anno <- read.table(paste0(.x, ".fuma.aux.txt"), header = T, sep="\t")
#   fuma_p <- fuma_res$P[match(target_cell_type, fuma_anno$cell_type)]
#   fuma_fdr <- fuma_res$FDR[match(target_cell_type, fuma_anno$cell_type)]
#   
#   
#   magma_res <-  read.table(paste0(.x,".magma.gsa.out"), header = T) %>% mutate(FDR = p.adjust(P, method="fdr"))
#   magma_anno <- read.table(paste0(.x,".magma.aux.txt"), header = T, sep="\t")
#   magma_p <- magma_res$P[match(target_cell_type, magma_anno$cell_type)]
#   magma_fdr <- magma_res$FDR[match(target_cell_type, magma_anno$cell_type)]
#   
#   seismic_res <- read.table(paste0(.x,".results.txt"), header = T, sep="\t")
#   seismic_p <- seismic_res$pvalue[match(target_cell_type, seismic_res$cell_type)]
#   seismic_fdr <- seismic_res$FDR[match(target_cell_type, seismic_res$cell_type)]
#   
#   # gene_anno <- read.table(.y, header = T)
#   # inf_genes <- read.table(paste0(.x,".inf_genes.txt"), header = T, sep="\t") %>%
#   #    left_join(gene_anno, by=c("gene" = "GENE")) %>% 
#   #    filter(is_influential)
#   # seismic_causal_ratio <- (nrow(filter(inf_genes, is_causal))/ nrow(inf_genes))
#   
#   # gene_anno <- read.table(.y,header = T)
#   # inf_genes <- read.table(paste0(.x,".inf_genes.txt"), header = T, sep="\t") %>%
#   #   left_join(gene_anno, by=c("gene" = "GENE")) %>% 
#   #   filter(is_influential)
#   # seismic_causal_ratio <- (nrow(filter(inf_genes, is_causal))/ nrow(inf_genes))
#   
#   if (file.exists(paste0(.x,".scdrs_res.csv"))){
#     scdrs_res <- read.table(paste0(.x,".scdrs_res.csv"), header = T, sep="\t") %>%
#       mutate(P = pnorm(assoc_mcz, lower.tail = F)) %>%
#       mutate(FDR = p.adjust(P, method="fdr"))
#     scdrs_p <-  scdrs_res$P[match(target_cell_type, scdrs_res$group)]
#     scdrs_fdr <-  scdrs_res$FDR[match(target_cell_type, scdrs_res$group)]
#   }else{
#     scdrs_p <- rep(NA,3)
#     scdrs_fdr <- rep(NA,3)
#   }
#   c(seismic_fdr, seismic_p, scdrs_fdr, scdrs_p, fuma_p, fuma_fdr, magma_p, magma_fdr)
# })%>% list_transpose() %>%
#   set_names(c("seismic_fdr.cell_1","seismic_fdr.cell_2","seismic_fdr.cell_3", "seismic_p.cell_1", "seismic_p.cell_2","seismic_p.cell_3", 
#               "scdrs_fdr.cell_1", "scdrs_fdr.cell_2", "scdrs_fdr.cell_3","scdrs_p.cell_1", "scdrs_p.cell_2","scdrs_p.cell_3",
#               "fuma_p.cell_1", "fuma_p.cell_2","fuma_p.cell_3","fuma_fdr.cell_1", "fuma_fdr.cell_2","fuma_fdr.cell_3",
#               "magma_p.cell_1","magma_p.cell_2","magma_p.cell_3", "magma_fdr.cell_1", "magma_fdr.cell_2", "magma_fdr.cell_3")) %>%
#   as_tibble() %>%
#   mutate(effect_size = strsplit(extended_para_df_multi$es_name, split="_") %>% map(~.x[2]) %>% unlist, 
#          scheme = c(rep("independent genes", nrow(extended_para_df_random)), rep("exclusive genes", nrow(extended_para_df_no_overlap)), rep("shared causal genes", nrow(extended_para_df_causal_overlap)), rep("shared other genes", nrow(extended_para_df_ct_overlap))),
#          gene_set = extended_para_df_multi$gs_name, expr_set = extended_para_df_multi$expr_name) %>%
#   select(contains("fdr"), effect_size, scheme, gene_set, expr_set) %>%
#   pivot_longer(cols = !effect_size & !scheme & !gene_set & !expr_set, names_prefix = "*_fdr\\.", names_to = "method.cell_type", values_to = "fdr") %>%
#   mutate(method = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[1]) %>% unlist, cell_type = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
#   mutate(method = recode(method, seismic_fdr = "seismic", scdrs_fdr = "scDRS", "magma_fdr" = "S-MAGMA", "fuma_fdr" = "FUMA"))
# 
# extended_power_df_multi <- extended_res_df_multi  %>% 
#   group_by(effect_size, gene_set, method, scheme) %>% 
#   summarise(power = length(which(fdr<=0.05))/n()) %>%
#   ungroup() %>%
#   distinct(effect_size, gene_set, method, scheme, power)%>%
#   group_by(effect_size, method, scheme) %>%
#   mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
#   distinct(effect_size, method, scheme, mean_power, sd_power)
#   
# 
# ggplot(extended_power_df_multi, aes(x = effect_size, y = mean_power,fill=method)) +
#   geom_bar(stat = "identity", position=position_dodge()) +
#   facet_grid(cols = vars(scheme)) +
#   scale_fill_manual(values = color_mapping_vec) +
#   theme_classic() +
#   ggtitle("Power of methods when multiple cell types are perturbed") + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(extended_power_df_multi, aes(x = effect_size, y = mean_power, color = method, group = method)) +
#   geom_ribbon(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power), fill = method), alpha = 0.3, color = NA) + 
#   geom_line(aes(y = mean_power)) +
#   facet_grid(cols = vars(scheme)) +
#   scale_color_manual(values = color_mapping_vec) +
#   scale_fill_manual(values = color_mapping_vec) +
#   theme_classic() +
#   ggtitle("Random disease gene sets sampled with Z-score weight") +
#   theme(plot.title = element_text(hjust = 0.5))


###### plot for multiple cell types ###### 500 genes #####
extended_para_df_file_random_500 <- here("data","expr","causal_sim","sc3_multi_ct_500","random","perturbed_expr", "parameter_df.txt")
extended_para_df_file_no_overlap_500 <- here("data","expr","causal_sim","sc3_multi_ct_500","no_overlap","perturbed_expr", "parameter_df.txt")
extended_para_df_file_causal_overlap_500 <- here("data","expr","causal_sim","sc3_multi_ct_500","causal_overlap","perturbed_expr", "parameter_df.txt")
extended_para_df_file_ct_overlap_500 <- here("data","expr","causal_sim","sc3_multi_ct_500","ct_overlap","perturbed_expr", "parameter_df.txt")

extended_para_df_random_500 <- read.table(extended_para_df_file_random_500,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
extended_para_df_no_overlap_500 <- read.table(extended_para_df_file_no_overlap_500,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
extended_para_df_causal_overlap_500 <- read.table(extended_para_df_file_causal_overlap_500,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
extended_para_df_ct_overlap_500 <- read.table(extended_para_df_file_ct_overlap_500,header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)

#extended_para_df <- rbind(extended_para_df_500, extended_para_df_1000, extended_para_df_2000, extended_para_df_5000)
extended_para_df_multi_500 <- rbind(extended_para_df_random_500, extended_para_df_no_overlap_500, extended_para_df_causal_overlap_500, extended_para_df_ct_overlap_500)

pr_curve_data_multi <- list()
extended_res_df_multi_500 <- map2(extended_para_df_multi_500$output_header, extended_para_df_multi_500$gene_anno_file, ~{
    target_cell_type <- c("target_cell_type_1", "target_cell_type_2", "target_cell_type_3")
    
    fuma_res <- read.table(paste0(.x, ".fuma.gsa.out"), header = TRUE) %>%
      mutate(FDR = p.adjust(P, method = "fdr"))
    fuma_anno <- read.table(paste0(.x, ".fuma.aux.txt"), header = TRUE, sep = "\t")
    fuma_ct_variable <- fuma_anno$encoded_name[match(target_cell_type, fuma_anno$cell_type)]
    fuma_p <- fuma_res$P[match(fuma_ct_variable, fuma_res$VARIABLE)]
    fuma_fdr <- fuma_res$FDR[match(fuma_ct_variable, fuma_res$VARIABLE)]
    fuma_ct_rank <- fuma_res %>%
      arrange(P) %>%
      mutate(rank = 1:n()) %>%
      filter(VARIABLE %in% fuma_ct_variable) %>%
      arrange(factor(VARIABLE, levels = fuma_ct_variable)) %>%
      pull(rank)
    
    magma_res <- read.table(paste0(.x, ".magma.gsa.out"), header = TRUE) %>%
      mutate(FDR = p.adjust(P, method = "fdr"))
    magma_anno <- read.table(paste0(.x, ".magma.aux.txt"), header = TRUE, sep = "\t")
    magma_ct_variable <- magma_anno$encoded_name[match(target_cell_type, magma_anno$cell_type)]
    magma_p <- magma_res$P[match(magma_ct_variable, magma_res$VARIABLE)]
    magma_fdr <- magma_res$FDR[match(magma_ct_variable, magma_res$VARIABLE)]
    magma_ct_rank <- magma_res %>%
      arrange(P) %>%
      mutate(rank = 1:n()) %>%
      filter(VARIABLE %in% magma_ct_variable) %>%
      arrange(factor(VARIABLE, levels = magma_ct_variable)) %>%
      pull(rank)
    
    seismic_res <- read.table(paste0(.x, ".results.txt"), header = TRUE, sep = "\t")
    seismic_p <- seismic_res$pvalue[match(target_cell_type, seismic_res$cell_type)]
    seismic_fdr <- seismic_res$FDR[match(target_cell_type, seismic_res$cell_type)]
    seismic_ct_rank <- seismic_res %>%
      arrange(pvalue) %>%
      mutate(rank = 1:n()) %>%
      filter(cell_type %in% target_cell_type) %>%
      arrange(factor(cell_type, levels = target_cell_type)) %>%
      pull(rank)
    
    gene_anno <- read.table(.y, header = TRUE) 
    if ("is_causal" %in% colnames(gene_anno)) {
      gene_anno <- map(1:3, ~select(gene_anno, "GENE", "is_causal")) %>%
        map(~set_names(.x, c("gene", "is_causal")))
    }else{
      gene_anno <- map(1:3, ~select(gene_anno, "GENE", paste0("is_causal_", .x))) %>%
        map(~set_names(.x, c("gene", "is_causal")))
    }
    
    inf_genes <- paste0(.x, ".", target_cell_type, ".inf_genes.txt") %>%
      map(~read.table(.x, header = TRUE, sep="\t")) %>%
      map2(gene_anno, ~left_join(.x, .y, by="gene")) %>% 
      map(~drop_na(.x, is_causal, dfbetas))
    seismic_pr <- map(inf_genes, ~PRROC::pr.curve(.x$dfbetas[.x$is_causal], .x$dfbetas[!.x$is_causal], curve = TRUE))
    pr_curve_data_multi[[.x]] <<- map(seismic_pr, ~.x$curve)
    
    causal_gene_num <- inf_genes %>% map(~filter(.x, is_causal) %>% nrow())
    tot_gene_num <- inf_genes %>% map(~nrow(.x))
    bg_rate <- map2(causal_gene_num, tot_gene_num, ~.x/.y)
    
    if (file.exists(paste0(.x, ".scdrs_res.csv"))){
      scdrs_res <- read.table(paste0(.x, ".scdrs_res.csv"), header = TRUE, sep = "\t") %>%
        mutate(P = pnorm(assoc_mcz, lower.tail = FALSE)) %>%
        mutate(FDR = p.adjust(P, method = "fdr"))
      scdrs_p <- scdrs_res$P[match(target_cell_type, scdrs_res$group)]
      scdrs_fdr <- scdrs_res$FDR[match(target_cell_type, scdrs_res$group)]
      scdrs_ct_rank <- scdrs_res %>%
        arrange(P) %>%
        mutate(rank = 1:n()) %>%
        filter(group %in% target_cell_type) %>%
        arrange(factor(group, levels = target_cell_type)) %>%
        pull(rank)
    } else {
      scdrs_p <- rep(NA, 3)
      scdrs_fdr <- rep(NA, 3)
      scdrs_ct_rank <- rep(NA, 3)
    }
    
    c(seismic_fdr, seismic_p, seismic_ct_rank, bg_rate, scdrs_fdr, scdrs_p, scdrs_ct_rank, fuma_p, fuma_fdr, fuma_ct_rank, magma_p, magma_fdr, magma_ct_rank)
  }
) %>% list_transpose() %>%
  set_names(c("seismic_fdr.cell_1","seismic_fdr.cell_2","seismic_fdr.cell_3", 
              "seismic_p.cell_1", "seismic_p.cell_2","seismic_p.cell_3", 
              "seismic_ct_rank.cell_1", "seismic_ct_rank.cell_2", "seismic_ct_rank.cell_3", 
              "seismic_causal_rate.cell_1", "seismic_causal_rate.cell_2", "seismic_causal_rate.cell_3", 
              "scdrs_fdr.cell_1", "scdrs_fdr.cell_2", "scdrs_fdr.cell_3",
              "scdrs_p.cell_1", "scdrs_p.cell_2","scdrs_p.cell_3", 
              "scdrs_ct_rank.cell_1", "scdrs_ct_rank.cell_2", "scdrs_ct_rank.cell_3", 
              "fuma_p.cell_1", "fuma_p.cell_2","fuma_p.cell_3",
              "fuma_fdr.cell_1", "fuma_fdr.cell_2","fuma_fdr.cell_3", 
              "fuma_ct_rank.cell_1", "fuma_ct_rank.cell_2", "fuma_ct_rank.cell_3", 
              "magma_p.cell_1","magma_p.cell_2","magma_p.cell_3", 
              "magma_fdr.cell_1", "magma_fdr.cell_2", "magma_fdr.cell_3", 
              "magma_ct_rank.cell_1", "magma_ct_rank.cell_2", "magma_ct_rank.cell_3" )) %>%
  as_tibble() %>%
  mutate(effect_size = strsplit(extended_para_df_multi_500$es_name, split="_") %>% map(~.x[2]) %>% unlist, 
         scheme = c(rep("independent genes", nrow(extended_para_df_random_500)), rep("exclusive genes", nrow(extended_para_df_no_overlap_500)), rep("shared causal genes", nrow(extended_para_df_causal_overlap_500)), rep("shared other genes", nrow(extended_para_df_ct_overlap_500))),
         gene_set = extended_para_df_multi_500$gs_name, expr_set = extended_para_df_multi_500$expr_name)

extended_power_df_multi_500 <- extended_res_df_multi_500 %>%
  select(contains("fdr"), effect_size, scheme, gene_set, expr_set) %>%
  pivot_longer(cols = !effect_size & !scheme & !gene_set & !expr_set, names_prefix = "*_fdr\\.", names_to = "method.cell_type", values_to = "fdr") %>%
  mutate(method = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[1]) %>% unlist, cell_type = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(method = recode(method, "seismic_fdr" = "seismic", "scdrs_fdr" = "scDRS", "magma_fdr" = "S-MAGMA", "fuma_fdr" = "FUMA"))

extended_power_df_multi_500 <- extended_power_df_multi_500  %>% 
  group_by(effect_size, gene_set, method, scheme) %>% 
  summarise(power = length(which(fdr<=0.05))/n()) %>%
  ungroup() %>%
  distinct(effect_size, gene_set, method, scheme, power)%>%
  group_by(effect_size, method, scheme) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  distinct(effect_size, method, scheme, mean_power, sd_power) %>% 
  mutate(causal_gene_num = 500)

# extended_power_df_multi <- extended_power_df_multi %>%
#   mutate(causal_gene_num = 400)

# extended_power_df_multi_all <- rbind(extended_power_df_multi_500, extended_power_df_multi)

ggplot(extended_power_df_multi_500  %>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power,fill=method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,linewidth=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(scheme)) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  theme(plot.title = element_text(hjust = 0.5))

# ggplot(extended_power_df_multi_all, aes(x = effect_size, y = mean_power, color = method, group = method)) +
#   geom_ribbon(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power), fill = method), alpha = 0.3, color = NA) + 
#   geom_line(aes(y = mean_power)) +
#   facet_grid(cols = vars(scheme), rows = vars(causal_gene_num)) +
#   scale_color_manual(values = color_mapping_vec) +
#   scale_fill_manual(values = color_mapping_vec) +
#   theme_classic() +
#   ggtitle("Random disease gene sets sampled with Z-score weight") +
#   theme(plot.title = element_text(hjust = 0.5))

#plot rank
extended_rank_df_multi_500 <- extended_res_df_multi_500 %>%
  select(contains("_ct_rank"), effect_size, scheme, gene_set, expr_set) %>%
  pivot_longer(cols = !effect_size & !scheme & !gene_set & !expr_set, names_prefix = "*_ct_rank\\.", names_to = "method.cell_type", values_to = "rank") %>%
  mutate(method = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[1]) %>% unlist, cell_type = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(method = recode(method, seismic_ct_rank = "seismic", scdrs_ct_rank = "scDRS", magma_ct_rank = "S-MAGMA", fuma_ct_rank = "FUMA"))

extended_rank_1_df_500 <- extended_rank_df_multi_500  %>% 
  group_by(effect_size, gene_set, scheme,  expr_set, method) %>%
  mutate(contain_top_rank = ifelse(any(rank==1), T,F)) %>%
  distinct(effect_size, scheme, gene_set, expr_set, contain_top_rank) %>%
  group_by(effect_size, gene_set, scheme, method) %>% 
  summarise(power = length(which(contain_top_rank))/length(which(!is.na(contain_top_rank))))  %>%
  group_by(effect_size,scheme, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  ungroup() %>%
  distinct(effect_size, scheme, mean_power, sd_power, method)  

ggplot(extended_rank_1_df_500%>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(scheme)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  #ggtitle("Power of ranking the target cell type as top 1") +
  theme(plot.title = element_text(hjust = 0.5))


#top 5
extended_rank_5_df_500 <- extended_rank_df_multi_500  %>% 
  group_by(effect_size, gene_set, scheme,  expr_set, method) %>%
  mutate(contain_top_rank = ifelse(all(rank<=5), T,F)) %>%
  distinct(effect_size, scheme, gene_set, expr_set, contain_top_rank) %>%
  group_by(effect_size, gene_set, scheme, method) %>% 
  summarise(power = length(which(contain_top_rank))/length(which(!is.na(contain_top_rank))))  %>%
  group_by(effect_size,scheme, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  ungroup() %>%
  distinct(effect_size, scheme, mean_power, sd_power, method)  

#top 5 supplementary plot
ggplot(extended_rank_5_df_500%>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(scheme)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  #ggtitle("Power of ranking the target cell type as top 1") +
  theme(plot.title = element_text(hjust = 0.5))

#plot for influential gene anlysis
#PRC plot: get the condition
pr_curve_data_multi <- pr_curve_data_multi %>% 
  map(~map(.x, ~set_colnames(.x, c("recall", "precision","score")))) %>% 
  map(~map(.x, ~as_tibble(.x))) %>% 
  map2(list_transpose(as.list(extended_res_df_multi_500 %>% select(seismic_causal_rate.cell_1, seismic_causal_rate.cell_2, seismic_causal_rate.cell_3))), ~map2(.x, .y, ~mutate(.x, precision = precision/.y))) %>%
  map(~purrr::reduce(.x, ~rbind(.x, .y))) %>%
  map(~approx(x = .x[["recall"]], y = .x[["precision"]], xout = seq(0.01, 1,0.01)))
  

plot_pr_curve_data_multi <- pr_curve_data_multi %>% 
  map(~as_tibble(.x)) %>% 
  map(~set_colnames(.x, c("recall","precision_fold"))) %>%
  map2(extended_res_df_multi_500$scheme, ~mutate(.x, scheme = .y)) %>%
  map2(extended_res_df_multi_500$effect_size, ~mutate(.x, effect_size = .y)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  group_by(scheme, effect_size, recall) %>%
  summarise( mean_precision_fold = mean(precision_fold), sd_precision_fold = sd(precision_fold)) %>%
  ungroup() %>%
  distinct(scheme, effect_size, recall, mean_precision_fold, sd_precision_fold) %>%
  mutate(effect_size = as.numeric(effect_size)) 

ggplot(plot_pr_curve_data_multi %>% mutate(effect_size = as.character(effect_size)), aes(x = recall, y = log2(mean_precision_fold), color = effect_size, group = effect_size)) +
  geom_line() +
  facet_grid(cols = vars(scheme)) +
  ylab("precision over random") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")) +
  theme_classic() +
  scale_color_manual(values = rev(viridis::viridis(6, alpha=0.8))) +
  guides(color = guide_legend(title = "log2(expression fold change)"))


###### plot for 10 cell types ###### 500 genes #####
many_ct_para_df_file_random <- here("data","expr","causal_sim","sc3_10_ct","random","perturbed_expr", "parameter_df.txt")
many_ct_para_df_file_no_overlap <- here("data","expr","causal_sim","sc3_10_ct","no_overlap","perturbed_expr", "parameter_df.txt")
many_ct_para_df_file_causal_overlap <- here("data","expr","causal_sim","sc3_10_ct","causal_overlap","perturbed_expr", "parameter_df.txt")
many_ct_para_df_file_ct_overlap <- here("data","expr","causal_sim","sc3_10_ct","ct_overlap","perturbed_expr", "parameter_df.txt")

many_ct_para_df_random <- read.table(many_ct_para_df_file_random, header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
many_ct_para_df_no_overlap <- read.table(many_ct_para_df_file_no_overlap, header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
#many_ct_para_df_causal_overlap <- read.table(many_ct_para_df_file_causal_overlap, header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
many_ct_para_df_ct_overlap <- read.table(many_ct_para_df_file_ct_overlap, header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)

#many_ct_para_df <- rbind(many_ct_para_df_random, many_ct_para_df_no_overlap, many_ct_para_df_causal_overlap, many_ct_para_df_ct_overlap)
many_ct_para_df <- rbind(many_ct_para_df_random, many_ct_para_df_no_overlap, many_ct_para_df_ct_overlap)

many_ct_res_df <- map2(many_ct_para_df$output_header, 1:nrow(many_ct_para_df), ~{
  
  if(.y %% 100 == 0){
    print(.y)
  }
  
  target_cell_type <- paste0("target_cell_type_", 1:10)
  
  fuma_res <- read.table(paste0(.x, ".fuma.gsa.out"), header = TRUE) %>%
    mutate(FDR = p.adjust(P, method = "fdr"))
  fuma_anno <- read.table(paste0(.x, ".fuma.aux.txt"), header = TRUE, sep = "\t")
  fuma_ct_variable <- fuma_anno$encoded_name[match(target_cell_type, fuma_anno$cell_type)]
  fuma_p <- fuma_res$P[match(fuma_ct_variable, fuma_res$VARIABLE)]
  fuma_fdr <- fuma_res$FDR[match(fuma_ct_variable, fuma_res$VARIABLE)]
  fuma_ct_rank <- fuma_res %>%
    arrange(P) %>%
    mutate(rank = 1:n()) %>%
    filter(VARIABLE %in% fuma_ct_variable) %>%
    arrange(factor(VARIABLE, levels = fuma_ct_variable)) %>%
    pull(rank)
  
  magma_res <- read.table(paste0(.x, ".magma.gsa.out"), header = TRUE) %>%
    mutate(FDR = p.adjust(P, method = "fdr"))
  magma_anno <- read.table(paste0(.x, ".magma.aux.txt"), header = TRUE, sep = "\t")
  magma_ct_variable <- magma_anno$encoded_name[match(target_cell_type, magma_anno$cell_type)]
  magma_p <- magma_res$P[match(magma_ct_variable, magma_res$VARIABLE)]
  magma_fdr <- magma_res$FDR[match(magma_ct_variable, magma_res$VARIABLE)]
  magma_ct_rank <- magma_res %>%
    arrange(P) %>%
    mutate(rank = 1:n()) %>%
    filter(VARIABLE %in% magma_ct_variable) %>%
    arrange(factor(VARIABLE, levels = magma_ct_variable)) %>%
    pull(rank)
  
  seismic_res <- read.table(paste0(.x, ".results.txt"), header = TRUE, sep = "\t")
  seismic_p <- seismic_res$pvalue[match(target_cell_type, seismic_res$cell_type)]
  seismic_fdr <- seismic_res$FDR[match(target_cell_type, seismic_res$cell_type)]
  seismic_ct_rank <- seismic_res %>%
    arrange(pvalue) %>%
    mutate(rank = 1:n()) %>%
    filter(cell_type %in% target_cell_type) %>%
    arrange(factor(cell_type, levels = target_cell_type)) %>%
    pull(rank)
  
  if (file.exists(paste0(.x, ".scdrs_res.csv"))){
    scdrs_res <- read.table(paste0(.x, ".scdrs_res.csv"), header = TRUE, sep = "\t") %>%
      mutate(P = pnorm(assoc_mcz, lower.tail = FALSE)) %>%
      mutate(FDR = p.adjust(P, method = "fdr"))
    scdrs_p <- scdrs_res$P[match(target_cell_type, scdrs_res$group)]
    scdrs_fdr <- scdrs_res$FDR[match(target_cell_type, scdrs_res$group)]
    scdrs_ct_rank <- scdrs_res %>%
      arrange(P) %>%
      mutate(rank = 1:n()) %>%
      filter(group %in% target_cell_type) %>%
      arrange(factor(group, levels = target_cell_type)) %>%
      pull(rank)
  } else {
    scdrs_p <- rep(NA, 10)
    scdrs_fdr <- rep(NA, 10)
    scdrs_ct_rank <- rep(NA, 10)
  }
  
  c(seismic_fdr, seismic_p, seismic_ct_rank, scdrs_fdr, scdrs_p, scdrs_ct_rank, fuma_p, fuma_fdr, fuma_ct_rank, magma_p, magma_fdr, magma_ct_rank)
}) %>% list_transpose() %>%
  set_names(c(paste0("seismic_fdr.cell_", 1:10), paste0("seismic_p.cell_", 1:10), paste0("seismic_ct_rank.cell_", 1:10),
              paste0("scdrs_fdr.cell_", 1:10), paste0("scdrs_p.cell_", 1:10), paste0("scdrs_ct_rank.cell_", 1:10),
              paste0("fuma_fdr.cell_", 1:10), paste0("fuma_p.cell_", 1:10), paste0("fuma_ct_rank.cell_", 1:10),
              paste0("magma_fdr.cell_", 1:10), paste0("magma_p.cell_", 1:10), paste0("magma_ct_rank.cell_", 1:10) )) %>%
  as_tibble() %>%
  mutate(effect_size = strsplit(many_ct_para_df$es_name, split="_") %>% map(~.x[2]) %>% unlist, 
         #scheme = c(rep("random", nrow(many_ct_para_df_random)), rep("no_overlap", nrow(many_ct_para_df_no_overlap)), rep("causal_overlap", nrow(many_ct_para_df_causal_overlap)), rep("ct_overlap", nrow(many_ct_para_df_ct_overlap))),
         scheme = c(rep("random", nrow(many_ct_para_df_random)), rep("no_overlap", nrow(many_ct_para_df_no_overlap)), rep("ct_overlap", nrow(many_ct_para_df_ct_overlap))),
         gene_set = many_ct_para_df$gs_name, expr_set = many_ct_para_df$expr_name) 

power_df_10_ct <- many_ct_res_df %>%
  select(contains("fdr"), effect_size, scheme, gene_set, expr_set) %>%
  pivot_longer(cols = !effect_size & !scheme & !gene_set & !expr_set, names_prefix = "*_fdr\\.", names_to = "method.cell_type", values_to = "fdr") %>%
  mutate(method = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[1]) %>% unlist, cell_type = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(method = recode(method, seismic_fdr = "seismic", scdrs_fdr = "scDRS", "magma_fdr" = "S-MAGMA", "fuma_fdr" = "FUMA"))

power_df_10_ct <- power_df_10_ct  %>% 
  group_by(effect_size, gene_set, method, scheme) %>% 
  summarise(power = length(which(fdr<=0.05))/n()) %>%
  ungroup() %>%
  distinct(effect_size, gene_set, method, scheme, power)%>%
  group_by(effect_size, method, scheme) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  distinct(effect_size, method, scheme, mean_power, sd_power) %>% 
  mutate(causal_gene_num = 500)

ggplot(power_df_10_ct  %>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power,fill=method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(scheme)) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  theme(plot.title = element_text(hjust = 0.5))

#rank
rank_df_10_ct <- many_ct_res_df %>%
  select(contains("_ct_rank"), effect_size, scheme, gene_set, expr_set) %>%
  pivot_longer(cols = !effect_size & !scheme & !gene_set & !expr_set, names_prefix = "*_ct_rank\\.", names_to = "method.cell_type", values_to = "rank") %>%
  mutate(method = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[1]) %>% unlist, cell_type = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(method = recode(method, seismic_ct_rank = "seismic", scdrs_ct_rank = "scDRS", magma_ct_rank = "S-MAGMA", fuma_ct_rank = "FUMA"))

rank_1_df_10_ct <- rank_df_10_ct  %>% 
  group_by(effect_size, gene_set, scheme,  expr_set, method) %>%
  mutate(contain_top_rank = ifelse(any(rank==1), T,F)) %>%
  distinct(effect_size, scheme, gene_set, expr_set, contain_top_rank) %>%
  group_by(effect_size, gene_set, scheme, method) %>% 
  summarise(power = length(which(contain_top_rank))/length(which(!is.na(contain_top_rank))))  %>%
  group_by(effect_size,scheme, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  ungroup() %>%
  distinct(effect_size, scheme, mean_power, sd_power, method)  

#top 1 supplementary plot
ggplot(rank_1_df_10_ct%>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(scheme)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  #ggtitle("Power of ranking the target cell type as top 1") +
  theme(plot.title = element_text(hjust = 0.5))
