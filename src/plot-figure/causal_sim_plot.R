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

#color
color_mapping_vec <- c("scDRS" ="#46B8DAFF" , "seismic" = "#D43F3AFF" , "S-MAGMA"="#EEA236FF" , "FUMA"="#5CB85CFF")

###### plot for extended ######
standard_para_df_file<- here("data","expr","causal_sim","sc3_standard","perturbed_expr", "parameter_df.txt")

standard_para_df <- read.table(standard_para_df_file,header=T, sep = " ") 

#load data
pr_curve_data <- list()

standard_res_df <- map2(standard_para_df$output_header, standard_para_df$gene_anno_file, ~{
  #FUMA results
  fuma_res <- read.table(paste0(.x, ".fuma.gsa.out"), header = T) %>% 
    mutate(FDR = p.adjust(P, method="fdr"))
  fuma_anno <- read.table(paste0(.x, ".fuma.aux.txt"), header = T, sep="\t")
  fuma_ct_variable <- fuma_anno$encoded_name[fuma_anno$cell_type == "target_cell_type"]
  fuma_p <- fuma_res$P[which(fuma_res$VARIABLE == fuma_ct_variable)]
  fuma_fdr <- fuma_res$FDR[which(fuma_res$VARIABLE == fuma_ct_variable)]
  fuma_ct_rank <- fuma_res %>% arrange(P) %>% mutate(rank = 1:n()) %>% filter(VARIABLE == fuma_ct_variable) %>% pull(rank)
  
  #S-MAGMA results
  magma_res <-  read.table(paste0(.x,".magma.gsa.out"), header = T) %>% 
    mutate(FDR = p.adjust(P, method="fdr"))
  magma_anno <- read.table(paste0(.x,".magma.aux.txt"), header = T, sep="\t")
  magma_ct_variable <- magma_anno$encoded_name[magma_anno$cell_type == "target_cell_type"]
  magma_p <- magma_res$P[which(magma_res$VARIABLE == magma_ct_variable)]
  magma_fdr <- magma_res$FDR[which(magma_res$VARIABLE == magma_ct_variable)]
  magma_ct_rank <- magma_res %>% arrange(P) %>% mutate(rank = 1:n()) %>% filter(VARIABLE == magma_ct_variable) %>% pull(rank)
  
  #scdrs results 
  scdrs_res <- read.table(paste0(.x,".scdrs_res.csv"), header = T, sep="\t") %>%
    mutate(P = pnorm(assoc_mcz, lower.tail = F)) %>%
    mutate(FDR = p.adjust(P, method="fdr"))
  scdrs_p <-  scdrs_res$P[which(scdrs_res$group=="target_cell_type")]
  scdrs_fdr <-  scdrs_res$FDR[which(scdrs_res$group=="target_cell_type")]
  scdrs_ct_rank <- scdrs_res %>% arrange(P) %>% mutate(rank = 1:n()) %>% filter(group == "target_cell_type") %>% pull(rank)
  
  #seismic results
  seismic_res <- read.table(paste0(.x,".results.txt"), header = T, sep="\t")
  seismic_p <- seismic_res$pvalue[which(seismic_res$cell_type=="target_cell_type")]
  seismic_fdr <- seismic_res$FDR[which(seismic_res$cell_type=="target_cell_type")]
  seismic_ct_rank <- seismic_res %>% arrange(pvalue) %>% mutate(rank = 1:n()) %>% filter(cell_type == "target_cell_type") %>% pull(rank)
  
  #influential gene analysis - pr curve
  gene_anno <- read.table(.y, header = T)
  inf_genes <- read.table(paste0(.x,".inf_genes.txt"), header = T, sep="\t") %>%
    left_join(gene_anno, by=c("gene" = "GENE")) %>% 
    drop_na(is_causal, dfbetas)
  pos_dfbetas <- inf_genes$dfbetas[inf_genes$is_causal]
  neg_dfbetas <- inf_genes$dfbetas[!inf_genes$is_causal]
  seismic_pr <- PRROC::pr.curve(pos_dfbetas, neg_dfbetas, curve = T)
  pr_curve_data[[.x]] <<- seismic_pr$curve
  
  #gene annotation information
  causal_gene_num <- length(pos_dfbetas)
  tot_gene_num <- nrow(inf_genes) 

  c(seismic_fdr, seismic_p, seismic_ct_rank, causal_gene_num, tot_gene_num, 
    scdrs_fdr, scdrs_p,scdrs_ct_rank, fuma_p, fuma_fdr,fuma_ct_rank, magma_p, magma_fdr,magma_ct_rank)
}) %>% list_transpose() %>%
  set_names(c("seismic_fdr","seismic_p", "seismic_ct_rank", "n_causal_gene","n_tot_gene",
              "scdrs_fdr", "scdrs_p","scdrs_ct_rank", "fuma_p", "fuma_fdr","fuma_ct_rank", "magma_p", "magma_fdr","magma_ct_rank")) %>%
  as_tibble() %>%
  mutate(effect_size = strsplit(standard_para_df$es_name, split="_") %>% map(~.x[2]) %>% unlist) %>%
  mutate(causal_gene_ratio = strsplit(standard_para_df$gr_name, split="_") %>% map(~.x[2]) %>% unlist) %>%
  mutate(data_set = standard_para_df$expr_name, 
         gene_set = standard_para_df$gs_name, 
         gene_num = standard_para_df$gene_num) 

standard_power_df <- standard_res_df  %>% 
  group_by(effect_size, causal_gene_ratio, gene_set) %>%
  mutate(seismic_ratio = length(which(seismic_fdr<=0.05))/length(which(!is.na(seismic_fdr))), 
         fuma_ratio = length(which(fuma_fdr<=0.05))/length(which(!is.na(fuma_fdr))),
         magma_ratio = length(which(magma_fdr<=0.05))/length(which(!is.na(magma_fdr))),
         scdrs_ratio = length(which(scdrs_fdr<=0.05))/length(which(!is.na(scdrs_fdr)))) %>%
  distinct( causal_gene_ratio, effect_size, gene_set, seismic_ratio, fuma_ratio, magma_ratio, scdrs_ratio) %>%
  set_colnames(c("causal_gene_ratio","effect_size","gene_set",  "seismic","FUMA", "S-MAGMA","scDRS")) %>%
  pivot_longer(cols = !effect_size & ! causal_gene_ratio & !gene_set, names_to = "method", values_to = "power") %>%
  group_by(effect_size, causal_gene_ratio, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  select(-gene_set) %>%
  distinct() %>%
  mutate(causal_gene_num = paste0(1000*as.numeric(causal_gene_ratio)," causal genes")) %>%
  mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA")))

#power plot
#figure main 
ggplot(standard_power_df %>% filter(causal_gene_num == "500 causal genes") , 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width = .4, linewidth = .3, position = position_dodge(.9)) + 
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  theme(plot.title = element_text(hjust = 0.5))

#figure supplementary
ggplot(standard_power_df %>% filter(causal_gene_ratio %in% c(0.2,0.4,0.5,0.6,0.8)), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4, linewidth= .3, position = position_dodge(.9)) + 
  facet_grid(cols = vars(causal_gene_num)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  theme(plot.title = element_text(hjust = 0.5))

#plot for the ranking
#top 1
standard_rank_1_df <- standard_res_df  %>% 
  group_by(effect_size, causal_gene_ratio, gene_set) %>%
  mutate(seismic_ratio = length(which(seismic_ct_rank<=1))/length(which(!is.na(seismic_ct_rank))), 
         fuma_ratio = length(which(fuma_ct_rank<=1))/length(which(!is.na(fuma_ct_rank))),
         magma_ratio = length(which(magma_ct_rank<=1))/length(which(!is.na(magma_ct_rank))),
         scdrs_ratio = length(which(scdrs_ct_rank<=1))/length(which(!is.na(scdrs_ct_rank)))) %>%
  distinct( causal_gene_ratio, effect_size,gene_set, seismic_ratio, fuma_ratio, magma_ratio, scdrs_ratio) %>%
  set_colnames(c("causal_gene_ratio","effect_size","gene_set",  "seismic","FUMA", "S-MAGMA","scDRS")) %>%
  pivot_longer(cols = !effect_size & !causal_gene_ratio & !gene_set, names_to = "method", values_to = "power") %>%
  group_by(effect_size, causal_gene_ratio, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  select(-gene_set) %>%
  distinct()  %>%
  mutate(causal_gene_num = paste0(1000*as.numeric(causal_gene_ratio), " causal genes")) %>%
  mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA")))

#top 1 supplementary plot
ggplot(standard_rank_1_df, 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(causal_gene_num)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  theme(plot.title = element_text(hjust = 0.5))

#top 5
standard_rank_5_df <- standard_res_df  %>% 
  group_by(effect_size, causal_gene_ratio, gene_set) %>%
  mutate(seismic_ratio = length(which(seismic_ct_rank<=5))/length(which(!is.na(seismic_ct_rank))), 
         fuma_ratio = length(which(fuma_ct_rank<=5))/length(which(!is.na(fuma_ct_rank))),
         magma_ratio = length(which(magma_ct_rank<=5))/length(which(!is.na(magma_ct_rank))),
         scdrs_ratio = length(which(scdrs_ct_rank<=5))/length(which(!is.na(scdrs_ct_rank)))) %>%
  distinct(causal_gene_ratio, effect_size,gene_set, seismic_ratio, fuma_ratio, magma_ratio, scdrs_ratio) %>%
  set_colnames(c("causal_gene_ratio","effect_size","gene_set", "seismic","FUMA", "S-MAGMA","scDRS")) %>%
  pivot_longer(cols = !effect_size & ! causal_gene_ratio & !gene_set, names_to = "method", values_to = "power") %>%
  group_by(effect_size, causal_gene_ratio, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  select(-gene_set) %>%
  distinct()  %>%
  mutate(causal_gene_num = paste0(1000*as.numeric(causal_gene_ratio), " causal genes")) %>%
  mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA")))

#top 1 supplementary plot
ggplot(standard_rank_5_df,
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(causal_gene_num)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  theme(plot.title = element_text(hjust = 0.5))

### Precision recall plot
pr_curve_data <- pr_curve_data %>% 
  map(~set_colnames(.x, c("recall", "precision","score"))) %>% 
  map(~as_tibble(.x))

plot_pr_curve_data <- pr_curve_data %>% 
  map2(standard_res_df$n_causal_gene / standard_res_df$n_tot_gene, ~mutate(.x, precision = precision/.y)) %>%
  split(standard_res_df$effect_size[standard_res_df$causal_gene_ratio == 0.5]) %>%
  map(~map(.x, ~approx(x = .x[["recall"]], y = .x[["precision"]], xout = seq(0.01,1,0.01)))) %>% #line smoothing
  map(~map(.x, ~as_tibble(.x))) %>%
  map(~map(.x, ~set_colnames(.x, c("recall", "precision_fold")))) %>% 
  map(~purrr::reduce(.x, ~rbind(.x, .y))) %>%
  map(~group_by(.x, recall)) %>%
  map(~mutate(.x, mean_precision_fold = mean(precision_fold), sd_precision_fold = sd(precision_fold))) %>%
  map(~ungroup(.x)) %>%
  map(~distinct(.x, recall, mean_precision_fold, sd_precision_fold)) %>%
  map2(names(.), ~mutate(.x, effect_size = as.numeric(.y))) %>%
  purrr::reduce(~rbind(.x, .y ))

ggplot(plot_pr_curve_data %>% mutate(effect_size = as.character(effect_size)), aes(x = recall, y = log2(mean_precision_fold), color = effect_size, group = effect_size)) +
  geom_line() +
  ylab("precision over random") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")) +
  theme_classic() +
  scale_color_manual(values = rev(viridis::viridis(6, alpha=0.8))) +
  guides(color = guide_legend(title = "log2(expression fold change)"))


#### plot for multiple cell types ####
extended_para_df_file_random <- here("data","expr","causal_sim","sc3_multi_ct","random","perturbed_expr", "parameter_df.txt")
extended_para_df_file_no_overlap <- here("data","expr","causal_sim","sc3_multi_ct","no_overlap","perturbed_expr", "parameter_df.txt")
extended_para_df_file_causal_overlap <- here("data","expr","causal_sim","sc3_multi_ct","causal_overlap","perturbed_expr", "parameter_df.txt")
extended_para_df_file_ct_overlap <- here("data","expr","causal_sim","sc3_multi_ct","ct_overlap","perturbed_expr", "parameter_df.txt")

extended_para_df_random <- read.table(extended_para_df_file_random, header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
extended_para_df_no_overlap <- read.table(extended_para_df_file_no_overlap, header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
extended_para_df_causal_overlap <- read.table(extended_para_df_file_causal_overlap, header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)
extended_para_df_ct_overlap <- read.table(extended_para_df_file_ct_overlap, header=T, sep = " ") %>% as_tibble() %>% mutate(gene_num = 1000)

#combine df
extended_para_df_multi <- rbind(extended_para_df_random, extended_para_df_no_overlap, extended_para_df_causal_overlap, extended_para_df_ct_overlap)

#load results
pr_curve_data_multi <- list()
extended_res_df_multi <- map2(extended_para_df_multi$output_header, extended_para_df_multi$gene_anno_file, ~{
  target_cell_type <- c("target_cell_type_1", "target_cell_type_2", "target_cell_type_3")
  
  #import FUMA results
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
  
  #import S-MAGMA results
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
  
  #import scDRS results
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

  #import seismic results
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
  
  #influential gene table
  inf_genes <- paste0(.x, ".", target_cell_type, ".inf_genes.txt") %>%
    map(~read.table(.x, header = TRUE, sep="\t")) %>%
    map2(gene_anno, ~left_join(.x, .y, by="gene")) %>% 
    map(~drop_na(.x, is_causal, dfbetas))
  seismic_pr <- map(inf_genes, ~PRROC::pr.curve(.x$dfbetas[.x$is_causal], .x$dfbetas[!.x$is_causal], curve = TRUE))
  pr_curve_data_multi[[.x]] <<- map(seismic_pr, ~.x$curve)
  
  causal_gene_num <- inf_genes %>% map(~filter(.x, is_causal) %>% nrow())
  tot_gene_num <- inf_genes %>% map(~nrow(.x))
  bg_rate <- map2(causal_gene_num, tot_gene_num, ~.x/.y)
  
  
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
  mutate(effect_size = strsplit(extended_para_df_multi$es_name, split="_") %>% map(~.x[2]) %>% unlist, 
         scheme = c(rep("independent genes", nrow(extended_para_df_random)), rep("exclusive genes", nrow(extended_para_df_no_overlap)), 
                    rep("shared causal genes", nrow(extended_para_df_causal_overlap)), rep("shared other genes", nrow(extended_para_df_ct_overlap))),
         gene_set = extended_para_df_multi$gs_name, expr_set = extended_para_df_multi$expr_name)

extended_power_df_multi <- extended_res_df_multi %>%
  select(contains("fdr"), effect_size, scheme, gene_set, expr_set) %>%
  pivot_longer(cols = !effect_size & !scheme & !gene_set & !expr_set, names_prefix = "*_fdr\\.", names_to = "method.cell_type", values_to = "fdr") %>%
  mutate(method = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[1]) %>% unlist, cell_type = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(method = recode(method, "seismic_fdr" = "seismic", "scdrs_fdr" = "scDRS", "magma_fdr" = "S-MAGMA", "fuma_fdr" = "FUMA"))

extended_power_df_multi <- extended_power_df_multi  %>% 
  group_by(effect_size, gene_set, method, scheme) %>% 
  summarise(power = length(which(fdr<=0.05))/n()) %>%
  ungroup() %>%
  distinct(effect_size, gene_set, method, scheme, power)%>%
  group_by(effect_size, method, scheme) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  distinct(effect_size, method, scheme, mean_power, sd_power) 

ggplot(extended_power_df_multi  %>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power,fill=method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,linewidth=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(scheme)) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  theme(plot.title = element_text(hjust = 0.5))


#plot rank
extended_rank_df_multi <- extended_res_df_multi %>%
  select(contains("_ct_rank"), effect_size, scheme, gene_set, expr_set) %>%
  pivot_longer(cols = !effect_size & !scheme & !gene_set & !expr_set, names_prefix = "*_ct_rank\\.", names_to = "method.cell_type", values_to = "rank") %>%
  mutate(method = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[1]) %>% unlist, cell_type = strsplit(method.cell_type, split = ".", fixed = T) %>% map(~.x[2]) %>% unlist) %>%
  mutate(method = recode(method, seismic_ct_rank = "seismic", scdrs_ct_rank = "scDRS", magma_ct_rank = "S-MAGMA", fuma_ct_rank = "FUMA"))

extended_rank_1_df_ <- extended_rank_df_multi  %>% 
  group_by(effect_size, gene_set, scheme,  expr_set, method) %>%
  mutate(contain_top_rank = ifelse(any(rank==1), T,F)) %>%
  distinct(effect_size, scheme, gene_set, expr_set, contain_top_rank) %>%
  group_by(effect_size, gene_set, scheme, method) %>% 
  summarise(power = length(which(contain_top_rank))/length(which(!is.na(contain_top_rank))))  %>%
  group_by(effect_size,scheme, method) %>%
  mutate(mean_power = mean(power, na.rm = T), sd_power = sd(power, na.rm = T)) %>%
  ungroup() %>%
  distinct(effect_size, scheme, mean_power, sd_power, method)  

ggplot(extended_rank_1_df_%>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
       aes(x = effect_size, y = mean_power, fill = method, group = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_power - sd_power), ymax = pmin(1, mean_power + sd_power)), width=.4,size=.3, position=position_dodge(.9)) + 
  facet_grid(cols = vars(scheme)) +
  scale_color_manual(values = color_mapping_vec) +
  scale_fill_manual(values = color_mapping_vec) +
  theme_classic() +
  xlab("log2(expression fold change)") +
  ylab("power") +
  theme(plot.title = element_text(hjust = 0.5))


#top 5
extended_rank_5_df_ <- extended_rank_df_multi  %>% 
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
ggplot(extended_rank_5_df_%>% mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))), 
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

#plot for influential gene analysis
#PRC plot
pr_curve_data_multi <- pr_curve_data_multi %>% 
  map(~map(.x, ~set_colnames(.x, c("recall", "precision","score")))) %>% 
  map(~map(.x, ~as_tibble(.x))) %>% 
  map2(list_transpose(as.list(extended_res_df_multi %>% select(seismic_causal_rate.cell_1, seismic_causal_rate.cell_2, seismic_causal_rate.cell_3))), ~map2(.x, .y, ~mutate(.x, precision = precision/.y))) %>%
  map(~purrr::reduce(.x, ~rbind(.x, .y))) %>%
  map(~approx(x = .x[["recall"]], y = .x[["precision"]], xout = seq(0.01, 1,0.01))) #smoothing

#pr curve plot 
plot_pr_curve_data_multi <- pr_curve_data_multi %>% 
  map(~as_tibble(.x)) %>% 
  map(~set_colnames(.x, c("recall","precision_fold"))) %>%
  map2(extended_res_df_multi$scheme, ~mutate(.x, scheme = .y)) %>%
  map2(extended_res_df_multi$effect_size, ~mutate(.x, effect_size = .y)) %>%
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

