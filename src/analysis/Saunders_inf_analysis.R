#Analysis of Saunders et al data set - influential gene analysis

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

library("clusterProfiler")
library("org.Hs.eg.db")

##### 2 calculate dfbetas #####
### PD
asso_tbl_pd = get_ct_asso(brain_sce, trait_name="PD_Nalls_part", asso_model="linear")
top_10_pd_ct = asso_tbl_pd %>% arrange(Pvalue)  %>% filter(1:n()<=10 & FDR<=0.05) %>% pull(cell_type) %>% set_names(.)
top_10_pd_dfbetas = map(top_10_pd_ct, ~gene_inf_measure(brain_sce, gene_zscore_df = gwas_zscore,cell_type = .x,trait_name = "PD_Nalls_part")) %>% 
  map(~mutate(.x, gene_symbol = mapIds(org.Hs.eg.db, keys = hsa_entrez, column = "SYMBOL",keytype = "ENTREZID")))
top_10_pd_inf_genes = top_10_pd_dfbetas %>%
  map(~filter(.x, PD_Nalls_part_z_stat>0 & influential)) %>%
  map(~pull(.x,hsa_entrez))
top_10_pd_inf_genes_neg = top_10_pd_dfbetas %>%
  map(~filter(.x, PD_Nalls_part_z_stat<0 & influential)) %>%
  map(~pull(.x,hsa_entrez)) 
top_10_pd_go = top_10_pd_inf_genes %>% 
  map2(~enrichGO(gene = .x , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr"))
top_10_pd_go_neg = top_10_pd_inf_genes_neg %>% 
  keep(~length(.x)>0) %>%
  map(~enrichGO(gene = .x , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr"))


###plot 
jaccard_mtx_pd = tibble(cell_type = top_10_pd_ct, genes = top_10_pd_inf_genes ) %>%
  cross_join(.,set_colnames(.,c("cell_type2","genes2"))) %>%
  rowwise() %>%
  mutate(jaccard_idx = jaccrd_idx(genes, genes2)) %>%
  ungroup() %>% 
  dplyr::select(cell_type, cell_type2, jaccard_idx) %>% 
  mutate(cell_type = factor(cell_type,levels=top_10_pd_ct), cell_type2 = factor(cell_type2,levels=top_10_pd_ct)) 
ggplot(jaccard_mtx_pd, aes(x=cell_type, y=cell_type2,fill=jaccard_idx)) + 
  geom_tile() + 
  viridis::scale_fill_viridis() +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  ylab("cell_type")

top_10_pd_go_plot = top_10_pd_go %>% 
  set_names(top_10_pd_ct) %>% 
  keep(~nrow(as.data.frame(.x))>0) %>%
  map2(names(.), ~dotplot(.x, title = paste0(.y)))
cowplot::plot_grid(plotlist = top_10_pd_go_plot, nrow=2)

top_10_pd_go_neg_plot = top_10_pd_go_neg %>% 
  keep(~nrow(as.data.frame(.x))>0) %>%
  map2(names(.), ~dotplot(.x, title = paste0(.y)))
cowplot::plot_grid(plotlist = top_10_pd_go_neg_plot, nrow=2)

#plot for genes
all_genes_pd = top_10_pd_dfbetas %>% 
  map(~filter(.x, influential)) %>%
  map2(names(.), ~mutate(.x,cell_type = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% 
  mutate(z_stat_value = ifelse(PD_Nalls_part_z_stat >=0, "positive","negative")) %>%
  mutate(cell_type = factor(cell_type, levels= top_10_pd_ct)) %>%
  group_by(cell_type, z_stat_value) %>% 
  summarise(Count = n(), .groups = 'drop') %>%
  add_row(cell_type = .$cell_type, z_stat_value = "negative", Count = 0) %>%
  group_by(cell_type, z_stat_value) %>% 
  filter(Count==max(Count)) %>% 
  ungroup()


ggplot(all_genes_pd, aes(x=cell_type,y=Count, fill=z_stat_value)) +  
  geom_bar(stat = "identity",position = "dodge") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  scale_fill_viridis_d()

### tau
asso_tbl_tau = get_ct_asso(brain_sce, trait_name="tau", asso_model="linear")
top_10_tau_ct = asso_tbl_tau %>% arrange(Pvalue)  %>% filter(1:n()<=10 & FDR<=0.05) %>% pull(cell_type)  %>% set_names(.)
top_10_tau_dfbetas = map(top_10_tau_ct, ~gene_inf_measure(brain_sce, gene_zscore_df = gwas_zscore,cell_type = .x,trait_name = "tau")) %>%
  map(~mutate(.x, gene_symbol = mapIds(org.Hs.eg.db, keys = hsa_entrez, column = "SYMBOL",keytype = "ENTREZID")))
top_10_tau_inf_genes = top_10_tau_dfbetas %>%
  map(~filter(.x, tau_z_stat>0 & influential)) %>%
  map(~pull(.x,hsa_entrez))
top_10_tau_inf_genes_neg = top_10_tau_dfbetas %>%
  map(~filter(.x, tau_z_stat<0 & influential)) %>%
  map(~pull(.x,hsa_entrez)) 
top_10_tau_go = top_10_tau_inf_genes %>% 
  map(~enrichGO(gene = .x , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr"))
top_10_tau_go_neg = top_10_tau_inf_genes_neg %>% 
  keep(~length(.x)>0) %>%
  map(~enrichGO(gene = .x , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr"))

###plot for overlap
jaccard_mtx_tau = tibble(cell_type = top_10_tau_ct, genes = top_10_tau_inf_genes ) %>%
  cross_join(.,set_colnames(.,c("cell_type2","genes2"))) %>%
  rowwise() %>%
  mutate(jaccard_idx = jaccrd_idx(genes, genes2)) %>%
  ungroup() %>% 
  dplyr::select(cell_type, cell_type2, jaccard_idx) %>% 
  mutate(cell_type = factor(cell_type,levels=top_10_tau_ct), cell_type2 = factor(cell_type2,levels=top_10_tau_ct)) 
ggplot(jaccard_mtx_tau, aes(x=cell_type, y=cell_type2,fill=jaccard_idx)) + 
  geom_tile() + 
  viridis::scale_fill_viridis() +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  ylab("cell_type")

top_10_tau_go_plot = top_10_tau_go %>% 
  set_names(top_10_tau_ct) %>% 
  keep(~nrow(as.data.frame(.x))>0) %>%
  map2(names(.), ~dotplot(.x, title = paste0(.y)))
cowplot::plot_grid(plotlist = top_10_tau_go_plot, nrow=2)

#plot for genes
all_genes_tau = top_10_tau_dfbetas %>% 
  map(~filter(.x, influential)) %>%
  map2(names(.), ~mutate(.x,cell_type = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% 
  mutate(z_stat_value = ifelse(tau_z_stat >=0, "positive","negative")) %>%
  mutate(cell_type = factor(cell_type, levels= top_10_tau_ct)) %>%
  group_by(cell_type, z_stat_value) %>% 
  summarise(Count = n(), .groups = 'drop') %>%
  add_row(cell_type = .$cell_type, z_stat_value = "negative", Count = 0) %>%
  group_by(cell_type, z_stat_value) %>% 
  filter(Count==max(Count)) %>% 
  ungroup()

ggplot(all_genes_tau, aes(x=cell_type,y=Count, fill=z_stat_value)) +  
  geom_bar(stat = "identity",position = "dodge") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
  scale_fill_viridis_d()

##write out
asso_tbl_tau %>% 
  arrange(Pvalue) %>%
  write.csv(file=here("results","inf_analysis","tau","enriched_cell_types.csv"))
top_10_tau_dfbetas %>% 
  map(~arrange(.x, -dfbetas)) %>%
  map2(names(.), ~write.csv(.x,here("results","inf_analysis","tau",paste0("dfbetas_",gsub(pattern = "/",replacement = " or ", fixed=T, x=.y),".csv"))))

### Kunkle AD
asso_tbl_Kunkle_AD = get_ct_asso(brain_sce, trait_name="Kunkleetal_2019", asso_model="linear")
top_10_Kunkle_AD_ct = asso_tbl_Kunkle_AD %>% arrange(Pvalue)  %>% filter(1:n()<=10 & FDR<=0.05) %>% pull(cell_type)  %>% set_names(.)
top_10_Kunkle_AD_dfbetas = map(top_10_Kunkle_AD_ct, ~gene_inf_measure(brain_sce, gene_zscore_df = gwas_zscore,cell_type = .x,trait_name = "Kunkleetal_2019")) %>%
  map(~mutate(.x, gene_symbol = mapIds(org.Hs.eg.db, keys = hsa_entrez, column = "SYMBOL",keytype = "ENTREZID")))
top_10_Kunkle_AD_inf_genes = top_10_Kunkle_AD_dfbetas %>%
  map(~filter(.x, Kunkleetal_2019_z_stat>0 & influential)) %>%
  map(~pull(.x,hsa_entrez))
top_10_Kunkle_AD_inf_genes_neg = top_10_Kunkle_AD_dfbetas %>%
  map(~filter(.x, Kunkleetal_2019_z_stat<0 & influential)) %>%
  map(~pull(.x,hsa_entrez)) 
top_10_Kunkle_AD_go = top_10_Kunkle_AD_inf_genes %>% 
  map(~enrichGO(gene = .x , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr"))

###write out
asso_tbl_Kunkle_AD %>% 
  arrange(Pvalue) %>%
  write.csv(file=here("results","inf_analysis","Kunkle_AD","enriched_cell_types.csv"))
top_10_Kunkle_AD_dfbetas %>% 
  map(~arrange(.x, -dfbetas)) %>%
  map2(names(.), ~write.csv(.x,here("results","inf_analysis","Kunkle_AD",paste0("dfbetas_",.y,".csv"))))



###plot for overlap
jaccard_mtx_Kunkle_AD = tibble(cell_type = top_10_Kunkle_AD_ct, genes = top_10_Kunkle_AD_inf_genes ) %>%
  cross_join(.,set_colnames(.,c("cell_type2","genes2"))) %>%
  rowwise() %>%
  mutate(jaccard_idx = jaccrd_idx(genes, genes2)) %>%
  ungroup() %>% 
  dplyr::select(cell_type, cell_type2, jaccard_idx) %>% 
  mutate(cell_type = factor(cell_type,levels=top_10_Kunkle_AD_ct), cell_type2 = factor(cell_type2,levels=top_10_Kunkle_AD_ct)) 
ggplot(jaccard_mtx_Kunkle_AD, aes(x=cell_type, y=cell_type2,fill=jaccard_idx)) + 
  geom_tile() + 
  viridis::scale_fill_viridis() +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  ylab("cell_type")

top_10_Kunkle_AD_go_plot = top_10_Kunkle_AD_go %>% 
  set_names(top_10_Kunkle_AD_ct) %>% 
  keep(~nrow(as.data.frame(.x))>0) %>%
  map2(names(.), ~dotplot(.x, title = paste0(.y)))
cowplot::plot_grid(plotlist = top_10_Kunkle_AD_go_plot, nrow=2)

#plot for genes
all_genes_Kunkle_AD = top_10_Kunkle_AD_dfbetas %>% 
  map(~filter(.x, influential)) %>%
  map2(names(.), ~mutate(.x,cell_type = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% 
  mutate(z_stat_value = ifelse(Kunkleetal_2019_z_stat >=0, "positive","negative")) %>%
  mutate(cell_type = factor(cell_type, levels= top_10_Kunkle_AD_ct)) %>%
  group_by(cell_type, z_stat_value) %>% 
  summarise(Count = n(), .groups = 'drop') %>%
  add_row(cell_type = .$cell_type, z_stat_value = "negative", Count = 0) %>%
  group_by(cell_type, z_stat_value) %>% 
  filter(Count==max(Count)) %>% 
  ungroup()

ggplot(all_genes_Kunkle_AD, aes(x=cell_type,y=Count, fill=z_stat_value)) +  
  geom_bar(stat = "identity",position = "dodge") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
  scale_fill_viridis_d()



### Plot for PD_PGCALZ
asso_tbl_PGCALZ = get_ct_asso(brain_sce, trait_name="PGCALZ2_2021", asso_model="linear")
top_10_PGCALZ_ct = asso_tbl_PGCALZ %>% arrange(Pvalue)  %>% filter(1:n()<=10 & FDR<=0.05) %>% pull(cell_type)  %>% set_names(.)
top_10_PGCALZ_dfbetas = map(top_10_PGCALZ_ct, ~gene_inf_measure(brain_sce, gene_zscore_df = gwas_zscore,cell_type = .x,trait_name = "PGCALZ2_2021")) %>%
  map(~mutate(.x, gene_symbol = mapIds(org.Hs.eg.db, keys = hsa_entrez, column = "SYMBOL",keytype = "ENTREZID")))
top_10_PGCALZ_inf_genes = top_10_PGCALZ_dfbetas %>%
  map(~filter(.x, PGCALZ2_2021_z_stat>0 & influential)) %>%
  map(~pull(.x,hsa_entrez))
top_10_PGCALZ_inf_genes_neg = top_10_PGCALZ_dfbetas  %>%
  map(~filter(.x, PGCALZ2_2021_z_stat<0 & influential)) %>%
  map(~pull(.x,hsa_entrez))  #no genes
top_10_PGCALZ_go = top_10_PGCALZ_inf_genes %>% 
  map(~enrichGO(gene = .x , OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",readable = T,pAdjustMethod = "fdr"))


###plot for overlap
jaccard_mtx_PGCALZ = tibble(cell_type = top_10_PGCALZ_ct, genes = top_10_PGCALZ_inf_genes ) %>%
  cross_join(.,set_colnames(.,c("cell_type2","genes2"))) %>%
  rowwise() %>%
  mutate(jaccard_idx = jaccrd_idx(genes, genes2)) %>%
  ungroup() %>% 
  dplyr::select(cell_type, cell_type2, jaccard_idx) %>% 
  mutate(cell_type = factor(cell_type,levels=top_10_PGCALZ_ct), cell_type2 = factor(cell_type2,levels=top_10_PGCALZ_ct)) 
ggplot(jaccard_mtx_PGCALZ, aes(x=cell_type, y=cell_type2,fill=jaccard_idx)) + 
  geom_tile() + 
  viridis::scale_fill_viridis(limits=c(0,1)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  ylab("cell_type")

top_10_PGCALZ_go_plot = top_10_PGCALZ_go %>% 
  set_names(top_10_PGCALZ_ct) %>% 
  keep(~nrow(as.data.frame(.x))>0) %>%
  map2(names(.), ~dotplot(.x, title = paste0(.y)))
cowplot::plot_grid(plotlist = top_10_PGCALZ_go_plot, nrow=2)

#plot for genes
all_genes_Kunkle_AD = top_10_Kunkle_AD_dfbetas %>% 
  map(~filter(.x, influential)) %>%
  map2(names(.), ~mutate(.x,cell_type = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% 
  mutate(z_stat_value = ifelse(Kunkleetal_2019_z_stat >=0, "positive","negative")) %>%
  mutate(cell_type = factor(cell_type, levels= top_10_Kunkle_AD_ct)) %>%
  group_by(cell_type, z_stat_value) %>% 
  summarise(Count = n(), .groups = 'drop') %>%
  add_row(cell_type = .$cell_type, z_stat_value = "negative", Count = 0) %>%
  group_by(cell_type, z_stat_value) %>% 
  filter(Count==max(Count)) %>% 
  ungroup()

ggplot(all_genes_Kunkle_AD, aes(x=cell_type,y=Count, fill=z_stat_value)) +  
  geom_bar(stat = "identity",position = "dodge") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
  scale_fill_viridis_d()


##### real plot for the enrichment
#### 
top_10_tau_go_plot_list = top_10_tau_go %>%
  keep(~nrow(as.data.frame(.x))>0) %>% 
  map(~as.data.frame(.x)) %>% 
  map(~group_by(.x,geneID)) %>%
  map(~mutate(.x,full_ID = paste0(ID, collapse = ","))) %>%
  map(~filter(.x, pvalue==min(pvalue))) %>%
  map(~ungroup(.x)) %>%
  map(~arrange(.x,p.adjust)) %>%
  map(~dplyr::slice(.x, 1:10)) %>%
  map(~mutate(.x, Description = factor(Description, levels = rev(Description)))) %>%
  map(~mutate(.x, tot_gene = as.numeric(gsub("[0-9]*/",replacement = "", x=GeneRatio)))) %>%
  map(~mutate(.x, `-log10(FDR)` = -log10(p.adjust))) %>%
  map( ~ggplot(data = .x, aes(x=Description, y=`-log10(FDR)`,size=Count,color=`-log10(FDR)`)) + 
         geom_point() + 
         coord_flip() +
         #scale_color_viridis_c(name = "FDR",limits = c(1,3),breaks=c(1,2,3),labels=c("0.05","0.01","0.001"),na.value = viridis::viridis(4)[4]) +
         scale_color_viridis_c(name = "FDR",
                               values = c(0, 0.4, 1),
                               limits = c(-log10(0.05),7),
                               breaks=c(-log10(0.05),3,7),
                               labels=c("0.05","1e-3","1e-7"),
                               na.value = viridis::viridis(5)[5]) +
         theme_classic() + 
         theme(axis.text.y = element_text(color="black", size=10)) + 
         ylim(0,4)  +
         scale_size_continuous(range = c(1,10), limits = c(1,20)) +
         xlab("GO terms") +
         #ggtitle(paste0("Top GO terms for influential genes of ",.y)) + 
         theme(aspect.ratio=4/3,plot.title = element_text(hjust=0, size=12))) %>%
  map2(names(.), ~{title =cowplot::ggdraw() + cowplot::draw_label(paste0("Top GO terms for influential genes of ",.y),  x=0.5, y=0.5, hjust=0.5) ;
  cowplot::plot_grid(title, .x, ncol = 1 ,rel_heights=c(0.025, 1))})

title = cowplot::ggdraw() + cowplot::draw_label("Top CSF tau associated cell types", fontface='bold', x=0.5, y=0.5, hjust=0.5)
cowplot::plot_grid(title, cowplot::plot_grid(plotlist = top_10_tau_go_plot_list ,nrow=2),ncol = 1 ,rel_heights=c(0.035, 1))

###write out
top_10_tau_go %>%
  keep(~nrow(as.data.frame(.x))>0) %>% 
  map(~as.data.frame(.x)) %>% 
  map(~group_by(.x,geneID)) %>%
  map(~mutate(.x,full_ID = paste0(ID, collapse = ","))) %>%
  map(~filter(.x, pvalue==min(pvalue))) %>%
  map(~ungroup(.x)) %>%
  map(~arrange(.x,p.adjust)) %>%
  map2(names(.),~write.csv(.x, file = here("results","inf_analysis","tau",paste0("GO_",gsub(pattern = "/",replacement = " or ", fixed=T, x=.y),".csv"))))

#### For PD
top_10_pd_go_plot_list = top_10_pd_go %>%
  keep(~nrow(as.data.frame(.x))>0) %>% 
  map(~as.data.frame(.x)) %>% 
  map(~group_by(.x,geneID)) %>%
  map(~mutate(.x,full_ID = paste0(ID, collapse = ","))) %>%
  map(~filter(.x, pvalue==min(pvalue))) %>%
  map(~ungroup(.x)) %>%
  map(~arrange(.x,p.adjust)) %>%
  map(~dplyr::slice(.x, 1:10)) %>%
  map(~mutate(.x, Description = factor(Description, levels = rev(Description)))) %>%
  map(~mutate(.x, tot_gene = as.numeric(gsub("[0-9]*/",replacement = "", x=GeneRatio)))) %>%
  map(~mutate(.x, `-log10(FDR)` = -log10(p.adjust))) %>%
  map( ~ggplot(data = .x, aes(x=Description, y=`-log10(FDR)`,size=Count,color=`-log10(FDR)`)) + 
         geom_point() + 
         coord_flip() +
         #scale_color_viridis_c(name = "FDR",limits = c(1,3),breaks=c(1,2,3),labels=c("0.05","0.01","0.001"),na.value = viridis::viridis(4)[4]) +
         scale_color_viridis_c(name = "FDR",limits = c(1,8),breaks=c(1,3,7),labels=c("0.05","1e-4","1e-8"),na.value = viridis::viridis(4)[4]) +
         theme_minimal() + 
         theme(axis.text.y = element_text(color="black", size=10)) + 
         ylim(0,4)  +
         scale_size_continuous(range = c(1,10), limits = c(1,20)) +
         xlab("GO terms") +
         #ggtitle(paste0("Top GO terms for influential genes of ",.y)) + 
         theme(aspect.ratio=4/3,plot.title = element_text(hjust=0, size=12))) %>%
  map2(names(.), ~{title =cowplot::ggdraw() + cowplot::draw_label(paste0("Top GO terms for influential genes of ",.y),  x=0.5, y=0.5, hjust=0.5) ;
  cowplot::plot_grid(title, .x, ncol = 1 ,rel_heights=c(0.025, 1))})

title = cowplot::ggdraw() + cowplot::draw_label("Top PD associated cell types", fontface='bold', x=0.5, y=0.5, hjust=0.5)
cowplot::plot_grid(title, cowplot::plot_grid(plotlist = top_10_pd_go_plot_list ,nrow=2),ncol = 1 ,rel_heights=c(0.035, 1))

#### For AD
top_10_Kunkle_go_plot_list = top_10_Kunkle_AD_go %>%
  keep(~nrow(as.data.frame(.x))>0) %>% 
  map(~as.data.frame(.x)) %>% 
  map(~group_by(.x,geneID)) %>%
  map(~mutate(.x,full_ID = paste0(ID, collapse = ","))) %>%
  map(~filter(.x, pvalue==min(pvalue))) %>%
  map(~ungroup(.x)) %>%
  map(~arrange(.x,p.adjust)) %>%
  map(~dplyr::slice(.x, 1:10)) %>%
  map(~mutate(.x, Description = factor(Description, levels = rev(Description)))) %>%
  map(~mutate(.x, tot_gene = as.numeric(gsub("[0-9]*/",replacement = "", x=GeneRatio)))) %>%
  map(~mutate(.x, `-log10(FDR)` = -log10(p.adjust))) %>%
  map( ~ggplot(data = .x, aes(x=Description, y=`-log10(FDR)`,size=Count,color=`-log10(FDR)`)) + 
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
         ylim(0,10)  +
         scale_size_continuous(range = c(1,10), limits = c(1,20)) +
         xlab("GO terms") +
         #ggtitle(paste0("Top GO terms for influential genes of ",.y)) + 
         theme(aspect.ratio=4/3,plot.title = element_text(hjust=0, size=12))) %>%
  map2(names(.), ~{title =cowplot::ggdraw() + cowplot::draw_label(paste0("Top GO terms for influential genes of ",.y),  x=0.5, y=0.5, hjust=0.5) ;
  cowplot::plot_grid(title, .x, ncol = 1 ,rel_heights=c(0.025, 1))})

title = cowplot::ggdraw() + cowplot::draw_label("Top AD associated cell types (Clinical AD GWAS of Kunkle et al.)", fontface='bold', x=0.5, y=0.5, hjust=0.5)
cowplot::plot_grid(title, cowplot::plot_grid(plotlist = top_10_Kunkle_go_plot_list ,nrow=2),ncol = 1 ,rel_heights=c(0.035, 1))

###write out
top_10_Kunkle_AD_go %>%
  keep(~nrow(as.data.frame(.x))>0) %>% 
  map(~as.data.frame(.x)) %>% 
  map(~group_by(.x,geneID)) %>%
  map(~mutate(.x,full_ID = paste0(ID, collapse = ","))) %>%
  map(~filter(.x, pvalue==min(pvalue))) %>%
  map(~ungroup(.x)) %>%
  map(~arrange(.x,p.adjust)) %>%
  map2(names(.),~write.csv(.x, file = here("results","inf_analysis","Kunkle_AD",paste0("GO_",.y,".csv"))))


top_10_PGCALZ_go_plot_list = top_10_PGCALZ_go %>%
  keep(~nrow(as.data.frame(.x))>0) %>% 
  map(~as.data.frame(.x)) %>% 
  map(~arrange(.x,p.adjust)) %>%
  map(~slice(.x, 1:10)) %>%
  map(~mutate(.x, Description = factor(Description, levels = rev(Description)))) %>%
  map(~mutate(.x, tot_gene = as.numeric(gsub("[0-9]*/",replacement = "", x=GeneRatio)))) %>%
  map(~mutate(.x, `-log10(FDR)` = -log10(p.adjust))) %>%
  map( ~ggplot(data = .x, aes(x=Description, y=`-log10(FDR)`,size=Count,color=`-log10(FDR)`)) + 
         geom_point() + 
         coord_flip() +
         scale_color_viridis_c(name = "FDR",limits = c(1,8),breaks=c(1,4,8),labels=c("0.05","1e-4","1e-8"),na.value = viridis::viridis(4)[4]) +
         theme_minimal() + 
         theme(axis.text.y = element_text(color="black", size=10)) + 
         ylim(0,10)  +
         scale_size_continuous(range = c(1,10), limits = c(1,20)) +
         xlab("GO terms") +
         #ggtitle(paste0("Top GO terms for influential genes of ",.y)) + 
         theme(aspect.ratio=4/3,plot.title = element_text(hjust=0, size=12))) %>%
  map2(names(.), ~{title =cowplot::ggdraw() + cowplot::draw_label(paste0("Top GO terms for influential genes of ",.y),  x=0.5, y=0.5, hjust=0.5) ;
  cowplot::plot_grid(title, .x, ncol = 1 ,rel_heights=c(0.025, 1))})

title = cowplot::ggdraw() + cowplot::draw_label("Top AD associated cell types (AD GWAS with proxy cases of PGCALZ2)", fontface='bold', x=0.5, y=0.5, hjust=0.5)
cowplot::plot_grid(title, cowplot::plot_grid(plotlist = top_10_PGCALZ_go_plot_list ,nrow=2),ncol = 1 ,rel_heights=c(0.035, 1))


#### plot for scatter plot
scatter_plot_Kunkle_AD = top_10_Kunkle_AD_dfbetas %>% 
  map2(names(.),~plot_gene_inf(.x, "gene_symbol",num_top_gene_label = 20, label_repel = T) +
         ggtitle(paste0("Influential genes of ", .y)) + theme(plot.title = element_text(hjust = 0.5)))
title = cowplot::ggdraw() + cowplot::draw_label("Influential genes for top AD associated cell types (Kunkle et al 2019 GWAS)", fontface='bold', x=0.5, y=0.5, hjust=0.5)
cowplot::plot_grid(title, cowplot::plot_grid(plotlist = scatter_plot_Kunkle_AD,nrow=2),ncol = 1 ,rel_heights=c(0.035, 1))

scatter_plot_PGCALZ_AD = top_10_PGCALZ_dfbetas %>% 
  map2(names(.),~plot_gene_inf(.x, "gene_symbol",num_top_gene_label = 20, label_repel = T) +
         ggtitle(paste0("Influential genes of ", .y)) + theme(plot.title = element_text(hjust = 0.5)))
title = cowplot::ggdraw() + cowplot::draw_label("AD GWAS with proxy cases of PGCALZ2", fontface='bold', x=0.5, y=0.5, hjust=0.5)
cowplot::plot_grid(title, cowplot::plot_grid(plotlist = scatter_plot_PGCALZ_AD ,nrow=2),ncol = 1 ,rel_heights=c(0.035, 1))


scatter_plot_tau_AD = top_10_tau_dfbetas %>% 
  map2(names(.),~plot_gene_inf(.x, "gene_symbol",num_top_gene_label = 20, label_repel = T) +
         ggtitle(paste0("Influential genes of ", .y)) + theme(plot.title = element_text(hjust = 0.5)))
title = cowplot::ggdraw() + cowplot::draw_label("Influential genes for top CSF tau associated cell types", fontface='bold', x=0.5, y=0.5, hjust=0.5)
cowplot::plot_grid(title, cowplot::plot_grid(plotlist = scatter_plot_tau_AD,nrow=2),ncol = 1 ,rel_heights=c(0.035, 1))

#### specifically for the other two terms
interested_genes = c("TRHR","TRHDE","MAP1LC3B","VPS26A","BECN1","ATP6V0E2","TEX264","RNF5","BACE2","MME","TRHR") 
respirtation_genes = as.data.frame(top_10_tau_go[[1]])$geneID[2] %>% strsplit(split="/") %>% unlist
autophagy_genes = as.data.frame(top_10_tau_go[[10]])$geneID[2] %>% strsplit(split="/") %>% unlist
ec_all_top5_genes = top_10_tau_dfbetas[[1]] %>% filter(influential) %>% arrange(-dfbetas) %>% slice(1:5) %>% pull(gene_symbol)
ec_2_top5_genes = top_10_tau_dfbetas[[10]] %>% filter(influential) %>% arrange(-dfbetas) %>% slice(1:5) %>% pull(gene_symbol)
ec_all_genes = intersect(c(interested_genes,respirtation_genes,autophagy_genes,ec_all_top10_genes),top_10_tau_dfbetas[[1]] %>% filter(hsa_entrez %in% top_10_tau_inf_genes[[1]]) %>% pull(gene_symbol))
ec_2_genes = intersect(c(interested_genes,respirtation_genes,autophagy_genes,ec_2_top10_genes),top_10_tau_dfbetas[[10]] %>% filter(hsa_entrez %in% top_10_tau_inf_genes[[10]]) %>% pull(gene_symbol))
ec_all_unqiue_genes = setdiff(ec_all_genes, ec_2_genes)
ec_2_unqiue_genes = setdiff(ec_2_genes, ec_all_genes)
common_genes = intersect(ec_2_genes, ec_all_genes)

ec_plot_df = top_10_tau_dfbetas[[1]] %>%
  mutate(influential = ifelse(influential & tau_z_stat>0, TRUE, FALSE)) %>%
  mutate(common = ifelse(gene_symbol %in% common_genes, TRUE, FALSE ))%>% 
  mutate(text_label = ifelse(gene_symbol %in% ec_all_genes, gene_symbol, ""))

ggplot(ec_plot_df, aes(x=specificity_score, y=tau_z_stat, color=common, label = text_label)) +
  geom_point(data = ec_plot_df %>% filter(!influential),aes(x=specificity_score, y = tau_z_stat), color="grey",alpha = 0.3 ) +
  geom_point(data = ec_plot_df %>% filter(influential),aes(x=specificity_score, y = tau_z_stat), color="red",alpha = 0.8 ) +
  geom_point(data = ec_plot_df %>% filter(text_label != ""),
             aes(x=specificity_score, y = tau_z_stat),
             color="black", shape=21, fill=NA,stroke=0.5) +
  theme_classic()+
  ggrepel::geom_text_repel( force  = 10,
                            alpha=1,
                            segment.square = F,
                            segment.inflect =T,
                            size=2,
                            max.overlaps = Inf,
                            segment.size=0.2) +
  scale_color_manual(values=c("black","blue2")) 

ec2_plot_df = top_10_tau_dfbetas[[10]] %>%
  mutate(influential = ifelse(influential & tau_z_stat>0, TRUE, FALSE)) %>%
  mutate(common = ifelse(gene_symbol %in% common_genes, TRUE, FALSE ))%>% 
  mutate(text_label = ifelse(gene_symbol %in% ec_2_genes, gene_symbol, ""))

ggplot(ec2_plot_df, aes(x=specificity_score, y=tau_z_stat, color=common, label = text_label)) +
  geom_point(data = ec2_plot_df %>% filter(!influential),aes(x=specificity_score, y = tau_z_stat), color="grey",alpha = 0.3 ) +
  geom_point(data = ec2_plot_df %>% filter(influential),aes(x=specificity_score, y = tau_z_stat), color="red",alpha = 0.8 ) +
  geom_point(data = ec2_plot_df %>% filter(text_label != ""),
             aes(x=specificity_score, y = tau_z_stat),
             color="black", shape=21, fill=NA,stroke=0.5) +
  theme_classic()+
  ggrepel::geom_text_repel( force  = 10,
                            alpha=1,
                            segment.square = F,
                            segment.inflect =T,
                            size=2,
                            max.overlaps = Inf,
                            segment.size=0.2) +
  scale_color_manual(values=c("black","blue2")) 
