#generate plots for results of the Saunders data set
##### load packages and results from scdrs/ours/fuma/magma #####
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

res_saunders <- list(fine_cluster = here("results","Saunders","fine_cluster","seismic","new_fine_cluster_res.txt"),
                subclass = here("results","Saunders","subclass","seismic","new_subclass_res.txt"),
                region_subclass = here("results","Saunders","region_subclass","seismic","new_region_subclass_res.txt"), 
                region_class = here("results","Saunders","region_class","seismic","new_region_classs_res.txt"),
                region_cluster = here("results","Saunders","region_cluster","seismic","new_region_cluster_res.txt")) %>%
  map(~read.table(.x, header = T, sep = "\t")) %>%
  map(~as_tibble(.x)) %>%
  map(~mutate(.x, cell_type = gsub(pattern = "\"", replacement = "", x = cell_type)))

res_scdrs = list(fine_cluster=here("results","Saunders","fine_cluster","scDRS"),
                 subclass = here("results","Saunders","subclass","scDRS"),
                 region_subclass = here("results","Saunders","region_subclass","scDRS"), 
                 region_class =here("results","Saunders","region_class","scDRS"),
                 region_cluster = here("results","Saunders","region_cluster","scDRS")) %>% 
  map(~list.files(.x, pattern = "^.*scdrs_group.*$",full.names = T)) %>% 
  map(~set_names(.x, str_extract(string = .x, pattern = "(?<=/)[^/]+$") %>% gsub(pattern = ".scdrs_group.*$", replacement = "", x=.) )) %>%
  map(~map(.x,~read.table(.x, header=T, sep="\t") %>% as_tibble())) %>%
  map(~map(.x, ~mutate(.x, P_val = pnorm(assoc_mcz, lower.tail = F)))) %>%
  map(~map(.x, ~dplyr::select(.x, any_of(c("group","P_val"))))) %>%
  map(~map2(.x, names(.x), ~set_colnames(.x, c("cell_type",.y)))) %>%
  map(~purrr::reduce(.x, ~left_join(.x, .y, by="cell_type"))) %>% 
  map(~mutate(.x, cell_type = gsub(pattern = "\"", replacement = "", fixed = T, x= cell_type))) %>% 
  map2(res_saunders, ~mutate(.x, cell_type = factor(cell_type, levels = .$cell_type))) %>% 
  map(~arrange(.x, cell_type))

###load fuma results
fuma_anno = list(file_cluster = here("data","expr","Saunders","new_Saunders.fine_cluster.fuma.aux.txt"), 
                 subclass = here("data","expr","Saunders","new_Saunders.subclass.fuma.aux.txt"),
                 region_subclass = here("data","expr","Saunders","new_Saunders.region_subclass.fuma.aux.txt"),
                 region_class = here("data","expr","Saunders","new_Saunders.region_class.fuma.aux.txt"),
                 region_cluster = here("data","expr","Saunders","new_Saunders.region_cluster.fuma.aux.txt")) %>%
  map(~read.table(.x, header = T, sep="\t") %>% as_tibble()) 

res_fuma = list(fine_cluster = here("results","new_Saunders","fine_cluster","FUMA"), 
                subclass = here("results","new_Saunders","subclass","FUMA"),
                region_subclass = here("results","new_Saunders","region_subclass","FUMA"), 
                region_class = here("results","new_Saunders","region_class","FUMA"),
                region_cluster = here("results","new_Saunders","region_cluster","FUMA")) %>%
  map(~list.files(.x, pattern = "*.gsa.out",full.names = T)) %>%
  map(~set_names(.x, str_extract(string=.x, pattern = "(?<=/)[^/]+$") %>%gsub(pattern = "\\.[A-Za-z.]+$", replacement = "", x=.) )) %>% 
  map(~map(.x, ~read.table(.x, header=T) %>% as_tibble())) %>% 
  map(~map(.x, ~dplyr::select(.x, any_of(c("VARIABLE","P"))))) %>%
  map(~map2(.x, names(.x), ~set_colnames(.x, c("cell_type",.y)))) %>%
  map(~purrr::reduce(.x, ~left_join(.x,.y,by="cell_type"))) %>% 
  map2(fuma_anno, ~left_join(.x %>% dplyr::rename(encoded_name = cell_type) , .y , by = "encoded_name")) %>%
  map(~dplyr::select(.x, -encoded_name) %>% relocate(cell_type) ) %>%
  map2(res_saunders, ~mutate(.x, cell_type = factor(cell_type, levels = .$cell_type))) %>%
  map(~arrange(.x, cell_type))

#load magma results
magma_anno = list(file_cluster = here("data","expr","Saunders","new_Saunders.fine_cluster.magma.aux.txt"), 
                 subclass = here("data","expr","Saunders","new_Saunders.subclass.magma.aux.txt"),
                 region_subclass = here("data","expr","Saunders","new_Saunders.region_subclass.magma.aux.txt"),
                 region_class = here("data","expr","Saunders","new_Saunders.region_class.magma.aux.txt"),
                 region_cluster = here("data","expr","Saunders","new_Saunders.region_cluster.magma.aux.txt")) %>%
  map(~read.table(.x, header = T, sep="\t") %>% as_tibble()) 

res_magma = list(fine_cluster = here("results","new_Saunders","fine_cluster","S-MAGMA"), 
                 subclass = here("results","new_Saunders","subclass","S-MAGMA"),
                 region_subclass = here("results","new_Saunders","region_subclass","S-MAGMA"), 
                 region_class = here("results","new_Saunders","region_class","S-MAGMA"),
                 region_cluster = here("results","new_Saunders","region_cluster","S-MAGMA"))  %>%
  map(~list.files(.x, pattern = "*.gsa.out",full.names = T)) %>%
  map(~set_names(.x, str_extract(string=.x, pattern = "(?<=/)[^/]+$") %>%gsub(pattern = "\\.[A-Za-z.]+$", replacement = "", x=.) )) %>% 
  map(~map(.x, ~read.table(.x, header=T) %>% as_tibble())) %>% 
  map(~map(.x, ~dplyr::select(.x, any_of(c("VARIABLE","P"))))) %>%
  map(~map2(.x, names(.x), ~set_colnames(.x, c("cell_type",.y)))) %>%
  map(~purrr::reduce(.x, ~left_join(.x,.y,by="cell_type"))) %>% 
  map2(magma_anno, ~left_join(.x %>% dplyr::rename(encoded_name = cell_type) , .y , by = "encoded_name")) %>%
  map(~dplyr::select(.x, -encoded_name) %>% relocate(cell_type) ) %>%
  map2(res_saunders, ~mutate(.x, cell_type = factor(cell_type, levels = .$cell_type))) %>%
  map(~arrange(.x, cell_type))

##### plot figure 4 #####
data_all = list("seismic" = res_saunders, "scDRS" = res_scdrs, "FUMA" = res_fuma, "S-MAGMA" = res_magma) %>% 
  map(~map(.x, ~pivot_longer(.x, !cell_type, names_to = "trait", values_to = "Pvalue"))) %>% 
  map(~map2(.x,names(.x), ~mutate(.x, granularity =.y ))) %>%
  map(~purrr::reduce(.x, ~rbind(.x,.y))) %>%
  map2(names(.), ~mutate(.x, method = .y)) %>%
  purrr::reduce(~rbind(.x, .y )) %>%
  group_by(method, granularity, trait) %>%
  mutate(FDR = p.adjust(Pvalue, method="fdr")) %>%
  ungroup()

#important information 
vn_neurons = c("SN.Neurons_SNc","SN.Neurons_VTA (ventral VTA)","SN.Neurons_SNc/VTA","SN.Neurons_Neurofilament state","SN.Neurons_VTA (dorsal VTA)",
               "Dopaminergic","SN.Dopaminergic","SN.Neurons","SN.DA neurons")

hline_data = data.frame(y = c(-log10(0.05),-log10(0.01)), type = factor(c("dotted","dashed")))

pd_plot_df = data_all %>% 
  filter(trait=="PD") %>% 
  mutate(neuron_type = ifelse(cell_type %in% vn_neurons, "DA neurons", ifelse(grepl(pattern="Endothelial|Oligodendrocytes|Polydendrocytes|Microglia_Macrophage|Astrocytes|Choroid_plexus|Mural|FibroblastLike",x=cell_type), "non-neurons", "other neurons"))) %>%
  mutate(neuron_type =factor(neuron_type, levels= c("DA neurons","other neurons", "non-neurons"))) %>% 
  group_by(method, granularity) %>%
  arrange(Pvalue) %>% 
  mutate(rank = 1:n()) %>% 
  ungroup %>% 
  mutate(text_label = ifelse(FDR<=0.01 |rank<=1  , as.character(cell_type), "")) %>%
  ungroup %>% 
  mutate(granularity = factor(granularity, levels = c("subclass","region_class","region_subclass","region_cluster","fine_cluster"))) %>%
  arrange(text_label) %>%
  mutate(method = factor(as.character(method),levels= c("S-MAGMA","FUMA","scDRS","seismic")))

ggplot(pd_plot_df, aes(x=method, y=-log10(FDR),color = neuron_type,alpha=neuron_type, label = text_label, size=neuron_type)) + 
  geom_jitter(position = position_jitter(seed = 1, width = 0.08)) + 
  facet_wrap(~granularity, ncol=1,
             labeller = labeller(granularity = c("subclass" = "subclass (n=14)","tissue_subclass" = "tissue_subclass (n=86)","fine_cluster" = "fine_cluster (n=214)", "tissue_class" = "tissue_class (n=75)", "tissue_cluster" = "tissue_cluster (n=108)"))) +
  theme_classic() +
  scale_size_manual(values=c(2,1,1)) +
  scale_color_manual(values=c(ggsci::pal_npg()(2),"grey60")) +
  scale_alpha_manual(values=c(0.9,0.5,0.4)) + 
  geom_hline(data= hline_data,aes(yintercept=y,linetype=type),linewidth=0.7,alpha=0.5) + 
  coord_flip() + 
  ggrepel::geom_text_repel( data = pd_plot_df %>% filter(text_label != ""),
                            force = 100,
                            position = ggpp::position_jitternudge(nudge.from = "jittered",
                                                                  width=0.08, y=4.5 + pd_plot_df %>% filter(text_label != "") %>% pull(FDR) %>% log10(), seed=1),
                            
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

#write out tables
pd_pval = pd_plot_df %>% 
  dplyr::select(cell_type, Pvalue, method,granularity) %>% 
  split(.$granularity) %>%
  map(~dplyr::select(.x,-granularity)) %>% 
  map(~pivot_wider(.x, names_from = method, values_from = Pvalue))
pd_fdr = pd_plot_df %>% 
  dplyr::select(cell_type, FDR, method,granularity) %>% 
  split(.$granularity) %>%
  map(~dplyr::select(.x,-granularity)) %>% 
  map(~pivot_wider(.x, names_from = method, values_from = FDR))
pd_pval %>% map2(names(.), ~write.csv(.x, here("figure","Saunders",paste0("PD_pval_",.y,"_association.csv"))))
pd_fdr %>% map2(names(.), ~write.csv(.x, here("figure","Saunders",paste0("PD_fdr_",.y,"_association.csv"))))

##### plot figure 5 #####
library(cowplot)
# for Kunkle AD 
kunkle_seismic = data_all %>% 
  filter(trait == "Kunkleetal_2019", method == "seismic", granularity=="fine_cluster") %>% 
  arrange(Pvalue) %>% 
  dplyr::slice(1:10) %>%
  mutate(cell_type = recode(cell_type, "ENT.Neurons_Entopeduncular Nucleus (Ent), GABA/Glu co-releasing" = "ENT neurons ENT, GABA/Glu")) %>%
  mutate(cell_type = gsub(pattern = "Microglia_", replacement = "microglia ", x=cell_type)) %>%
  mutate(cell_type = gsub(pattern = "Macrophage", replacement = "macrophage", x=cell_type)) %>%
  mutate(cell_type = gsub(pattern = "Neurons", replacement = "neurons", x=cell_type)) %>%
  mutate(cell_type = gsub(pattern = ".", replacement = " ",fixed = T, x=cell_type)) %>%
  mutate(cell_type = factor(cell_type,levels=rev(cell_type))) %>% 
  mutate(FDR = -log10(FDR)) 

kunkle_label_sig = kunkle_seismic  %>% dplyr::filter(FDR>= - log10(0.05) )%>% mutate(FDR = FDR + 0.3)

kunkle_ct = ggplot(kunkle_seismic, aes(x=cell_type, y = FDR)) + 
  ggplot2::geom_bar(stat = "identity",fill = "#D43F3AFF",alpha=0.8,width = 0.5)+ 
  coord_flip() +
  theme_minimal() +
  ylim(0,3.5) + 
  ggplot2::geom_text(data = kunkle_label_sig, label = "*") +
  ylab("-log10(FDR)") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") 

title_kunkle = ggdraw() + draw_label("Top AD-associated cell types, clinical AD GWAS", fontface='bold', x=0.5, y=0.5, hjust=0.5)
plot_grid(title_kunkle,kunkle_seismic,ncol=1,rel_heights = c(0.1,1))

# for tau
tau_seismic = data_all %>% 
  filter(trait == "tau", method == "seismic", granularity=="fine_cluster") %>% 
  arrange(Pvalue) %>% 
  dplyr::slice(1:10) %>%
  mutate(cell_type = factor(cell_type,levels=rev(cell_type))) %>% 
  mutate(FDR = -log10(FDR)) 

tau_label_sig = tau_seismic  %>% dplyr::filter(FDR>= - log10(0.05) )%>% mutate(FDR = FDR + 0.3)

tau_ct = ggplot(tau_seismic, aes(x=cell_type, y = FDR)) + 
  ggplot2::geom_bar(stat = "identity",fill = "#D43F3AFF",alpha=0.8,width = 0.5)+ 
  coord_flip() +
  theme_minimal() +
  ylim(0,3.5) + 
  ggplot2::geom_text(data = tau_label_sig, label = "*") +
  ylab("-log10(FDR)") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") 

title_tau = ggdraw() + draw_label("Top tau level-associated cell types, CSF tau GWAS", fontface='bold', x=0.5, y=0.5, hjust=0.5)
plot_grid(title_tau, tau_ct, ncol=1, rel_heights = c(0.1,1))

##import the influential gene analysis data and plot
library(org.Hs.eg.db)
Kunkle_hc_micro_dfbetas = read.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_dfbetas.txt"), header = T, sep = "\t") %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, keys = as.character(hsa_entrez), column = "SYMBOL",keytype = "ENTREZID"))
Kunkle_hc.mg_go = read.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_go.txt"), header=T, sep="\t")

tau_hc.ec_dfbetas = read.table(here("results","Saunders","inf_analysis","Kunkle_hc_ec_dfbetas.txt"), header = T, sep="\t") %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, keys = as.character(hsa_entrez), column = "SYMBOL",keytype = "ENTREZID"))
tau_hc.ec2_dfbetas = read.table(here("results","Saunders","inf_analysis","Kunkle_hc_ec2_dfbetas.txt"), header = T, sep="\t") %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, keys = as.character(hsa_entrez), column = "SYMBOL",keytype = "ENTREZID"))
tau_hc.ec_go = read.table(here("results","Saunders","inf_analysis","Kunkle_hc_ec_go.txt"), header = T, sep="\t")
tau_hc.ec2_go = read.table(here("results","Saunders","inf_analysis","Kunkle_hc_ec2_go.txt"), header = T, sep="\t")

#plot gene distribution
seismicGWAS::plot_gene_inf(Kunkle_hc_micro_dfbetas, gene_label_column = "gene_symbol", num_top_gene_label = 20)
seismicGWAS::plot_gene_inf(tau_hc.ec_dfbetas, gene_label_column = "gene_symbol", 
                           extra_labels = c("COX11","DCUN1D4","MRPL13","BECN1","MAP1LC3B","VPS26A","NDUFS2","NDUFA4","TBR1","ARPC1A","SYT10","DIXDC1","FEZF2","NPPC","NXPH3"))
seismicGWAS::plot_gene_inf(tau_hc.ec2_dfbetas,gene_label_column = "gene_symbol", 
                           extra_labels = c("GRAMD4","BECN1","MAP1LC3B", "VPS26A", "NDUFS2", "TRHR", "NTSR1", "CBLN4", "DCN", "PXDN", "RSPO2", "TMEM255B") )

#plot gene ontology plot
#mg with AD gwas: top influential genes
Kunkle_hc.mg_go.plot = Kunkle_hc.mg_go %>%
  group_by(geneID) %>%
  mutate(full_ID = paste0(ID, collapse = ",")) %>%
  filter(pvalue == min(pvalue)) %>%
  ungroup %>% 
  arrange(p.adjust ) %>% 
  dplyr::slice(1:10) %>% 
  mutate(Description = factor(Description, levels = rev(Description))) %>%
  mutate(`-log10(FDR)` = -log10(p.adjust)) %>%
  ggplot( aes(x=Description, y=`-log10(FDR)`, size=Count, color = `-log10(FDR)`)) +
    geom_point() + 
    coord_flip() +
    scale_color_viridis_c(name = "FDR",
                        values = c(0, 0.4, 1),
                        limits = c(-log10(0.05),7),
                        breaks=c(-log10(0.05),3,7),
                        labels=c("0.05","1e-3","1e-7"),
                        na.value = viridis::viridis(5)[5]) +
    theme_classic() + 
    theme(axis.text.y = element_text(color="black", size=10)) + 
    scale_size_continuous(range = c(1,10), limits = c(1,20)) +
    xlab("GO terms") +
    theme(aspect.ratio=4/3,plot.title = element_text(hjust=0, size=12))

#plot gene ontology plot
#ec with tau 
tau_hc.ec_go.plot = tau_hc.ec_go %>%
  group_by(geneID) %>%
  mutate(full_ID = paste0(ID, collapse = ",")) %>%
  filter(pvalue == min(pvalue)) %>%
  ungroup %>% 
  arrange(p.adjust ) %>% 
  dplyr::slice(1:10) %>% 
  mutate(Description = factor(Description, levels = rev(Description))) %>%
  mutate(`-log10(FDR)` = -log10(p.adjust)) %>%
  ggplot( aes(x=Description, y=`-log10(FDR)`, size=Count, color = `-log10(FDR)`)) +
  geom_point() + 
  coord_flip() +
  scale_color_viridis_c(name = "FDR",
                        values = c(0, 0.4, 1),
                        limits = c(-log10(0.05),7),
                        breaks=c(-log10(0.05),3,7),
                        labels=c("0.05","1e-3","1e-7"),
                        na.value = viridis::viridis(5)[5]) +
  theme_classic() + 
  theme(axis.text.y = element_text(color="black", size=10)) + 
  scale_size_continuous(range = c(1,10), limits = c(1,20)) +
  xlab("GO terms") +
  theme(aspect.ratio=4/3,plot.title = element_text(hjust=0, size=12))

#plot gene ontology plot
#ec2 with tau
tau_hc.ec2_go.plot = tau_hc.ec2_go  %>%
  group_by(geneID) %>%
  mutate(full_ID = paste0(ID, collapse = ",")) %>%
  filter(pvalue == min(pvalue)) %>%
  ungroup %>% 
  arrange(p.adjust ) %>% 
  dplyr::slice(1:10) %>% 
  mutate(Description = factor(Description, levels = rev(Description))) %>%
  mutate(`-log10(FDR)` = -log10(p.adjust)) %>%
  ggplot( aes(x=Description, y=`-log10(FDR)`, size=Count, color = `-log10(FDR)`)) +
  geom_point() + 
  coord_flip() +
  scale_color_viridis_c(name = "FDR",
                        values = c(0, 0.4, 1),
                        limits = c(-log10(0.05),7),
                        breaks=c(-log10(0.05),3,7),
                        labels=c("0.05","1e-3","1e-7"),
                        na.value = viridis::viridis(5)[5]) +
  theme_classic() + 
  theme(axis.text.y = element_text(color="black", size=10)) + 
  scale_size_continuous(range = c(1,10), limits = c(1,20)) +
  xlab("GO terms") +
  theme(aspect.ratio=4/3,plot.title = element_text(hjust=0, size=12))
