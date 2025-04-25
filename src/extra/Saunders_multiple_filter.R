#Analysis of Saunders et al data set

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

#load data
load(here("data","expr","Saunders","Saunders_processed.rda") )

##### specificity score at different cell number filter #####
all_sscore <- c(10, 20, 50) %>%
  set_names(c("10", "20", "50")) %>%
  map(~{ct_filter = .x; map(c("fine_cluster", "subclass", "region_class","region_subclass", "region_cluster"), ~calc_specificity(brain_sce, ct_label_col = .x, min_ct_size = ct_filter, min_avg_exp_ct = 0.01))}) %>%
  map(~map(.x, ~translate_gene_ids(.x, from = "mmu_symbol")))


#enrichment
gwas_zscore_pd <- here("data","gwas","brain_gwas","zscore","PD.genes.out")

all_results <- all_sscore %>%
  map(~map(.x, ~get_ct_trait_associations(sscore = .x, magma = gwas_zscore_pd)))


##### plot results ####
vn_neurons <- c("SN.Neurons_SNc","SN.Neurons_VTA (ventral VTA)","SN.Neurons_SNc/VTA","SN.Neurons_Neurofilament state","SN.Neurons_VTA (dorsal VTA)",
                "Dopaminergic","SN.Dopaminergic","SN.Neurons","SN.DA neurons")

hline_data <- data.frame(y = c(-log10(0.05),-log10(0.01)), type = factor(c("dotted","dashed")))

all_results <- all_results %>%
  map(~set_names(.x, c("fine_cluster", "subclass", "region_class","region_subclass", "region_cluster"))) %>%
  map(~map2(.x, names(.x), ~mutate(.x, granularity = .y))) %>%
  map(~purrr::reduce(.x, ~rbind(.x, .y))) %>%
  map2(names(.), ~mutate(.x, cell_number_filter = .y)) %>%
  purrr::reduce(~rbind(.x, .y))

all_results <- all_results %>%
  mutate(neuron_type = ifelse(cell_type %in% vn_neurons, "DA neurons", ifelse(grepl(pattern="Endothelial|Oligodendrocytes|Polydendrocytes|Microglia_Macrophage|Astrocytes|Choroid_plexus|Mural|FibroblastLike",x=cell_type), "non-neurons", "other neurons"))) %>%
  mutate(neuron_type =factor(neuron_type, levels= c("DA neurons","other neurons", "non-neurons"))) %>%
  group_by(cell_number_filter, granularity) %>%
  arrange(pvalue) %>% 
  mutate(rank = 1:n()) %>% 
  ungroup %>% 
  mutate(text_label = ifelse(FDR<=0.01 |rank<=1  , as.character(cell_type), "")) %>%
  ungroup %>% 
  mutate(granularity = factor(granularity, levels = c("subclass","region_class","region_subclass","region_cluster","fine_cluster"))) %>%
  arrange(text_label)

ggplot(all_results, aes(x=cell_number_filter, y=-log10(FDR), color = neuron_type, alpha=neuron_type, label = text_label, size=neuron_type)) + 
  geom_jitter(position = position_jitter(seed = 1, width = 0.08)) + 
  facet_wrap(~granularity, ncol=1,
             labeller = labeller(granularity = c("subclass" = "subclass (n=14)","tissue_subclass" = "tissue_subclass (n=86)","fine_cluster" = "fine_cluster (n=214)", "tissue_class" = "tissue_class (n=75)", "tissue_cluster" = "tissue_cluster (n=108)"))) +
  theme_classic() +
  scale_size_manual(values=c(2,1,1)) +
  scale_color_manual(values=c(ggsci::pal_npg()(2),"grey60")) +
  scale_alpha_manual(values=c(0.9,0.5,0.4)) + 
  geom_hline(data= hline_data,aes(yintercept=y,linetype=type),linewidth=0.7,alpha=0.5) + 
  coord_flip() + 
  ggrepel::geom_text_repel( data = all_results %>% filter(text_label != ""),
                            force = 100,
                            position = ggpp::position_jitternudge(nudge.from = "jittered",
                                                                  width=0.08, y=4.5 + all_results %>% filter(text_label != "") %>% pull(FDR) %>% log10(), seed=1),
                            
                            color="black",alpha=1,
                            segment.square = F,
                            segment.inflect =T,
                            size=2,
                            max.overlaps = Inf,
                            box.padding = 0.5,
                            segment.size = 0.2) +
  scale_linetype_manual(values = c("dashed","dotted"), 
                        labels = c("FDR<0.01", "FDR<0.05"),
                        name = "significance") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank())
