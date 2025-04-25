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

#load results
Kunkle_hc.micro_dfbetas <- read.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_dfbetas.txt"), header = T, sep = "\t") %>%
  as_tibble() %>%
  set_colnames(c("Entrez", "symbol", "specificity", "zscore", "dfbetas", "is_influential"))
Kunkle_hc.mg_go <- read.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_go.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(FDR)) %>%
  head(n=10) %>% 
  mutate(GO.name = factor(GO.name, levels = rev(GO.name)))

tau_hc.ec_dfbetas <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec_dfbetas.txt"), header = T, sep = "\t") %>%
  as_tibble()%>%
  set_colnames(c("Entrez", "symbol", "specificity", "zscore", "dfbetas", "is_influential"))
tau_hc.ec_go <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec_go.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(FDR)) %>% 
  head(n=10) %>%
  mutate(GO.name = factor(GO.name, levels = rev(GO.name)))
tau_hc.ec_kegg <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec_kegg.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(FDR)) %>% 
  head(n=10) %>%
  mutate(KEGG.name = factor(KEGG.name, levels = rev(KEGG.name)))

tau_hc.ec2_dfbetas <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec2_dfbetas.txt"), header = T, sep = "\t") %>%
  as_tibble() %>%
  set_colnames(c("Entrez", "symbol", "specificity", "zscore", "dfbetas", "is_influential"))
tau_hc.ec2_go <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec2_go.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(FDR)) %>% 
  head(n=10) %>%
  mutate(GO.name = factor(GO.name, levels = rev(GO.name)))
tau_hc.ec2_kegg <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec2_kegg.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(FDR)) %>% 
  head(n=10) %>%
  mutate(KEGG.name = factor(KEGG.name, levels = rev(KEGG.name)))

tau_hc.ca1_dfbetas <- read.table(here("results","Saunders","inf_analysis","tau_hc_ca1_dfbetas.txt"), header = T, sep = "\t") %>%
  as_tibble() %>%
  set_colnames(c("Entrez", "symbol", "specificity", "zscore", "dfbetas", "is_influential"))
tau_hc.ca1_go <- read.table(here("results","Saunders","inf_analysis","tau_hc_ca1_go.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(FDR)) %>% 
  head(n=10) %>%
  mutate(GO.name = factor(GO.name, levels = rev(GO.name)))
tau_hc.ca1_kegg <- read.table(here("results","Saunders","inf_analysis","tau_hc_ca1_kegg.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(FDR)) %>% 
  head(n=10) %>%
  mutate(KEGG.name = factor(KEGG.name, levels = rev(KEGG.name)))

#plot for influential genes
seismicGWAS::plot_inf_genes(Kunkle_hc.micro_dfbetas , gene_col = "genename")

seismicGWAS::plot_inf_genes(tau_hc.ec2_dfbetas, gene_col = "genename")

seismicGWAS::plot_inf_genes(tau_hc.ec_dfbetas, gene_col = "genename")

#venn diagram
library("ggvenn")
tau_hc.ec_inf_genes <- tau_hc.ec_dfbetas %>% filter(is_influential, zscore>0)  %>% pull(Entrez)
tau_hc.ec2_inf_genes <- tau_hc.ec2_dfbetas %>% filter(is_influential, zscore>0) %>% pull(Entrez)
tau_hc.ca1_inf_genes <- tau_hc.ca1_dfbetas %>% filter(is_influential, zscore>0) %>% pull(Entrez)

ggvenn(list(`CA1` = tau_hc.ca1_inf_genes, `entorhinal cortex` = tau_hc.ec_inf_genes, `entorhinal cortex layer II` = tau_hc.ec2_inf_genes), 
       text_size = 5,set_name_size=5,show_percentage = F,fill_color = c("#46B8DAFF",  "#D43F3AFF", "#EEA236FF"), fill_alpha = 0.4) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4))


#plot for GO enrichment
ggplot(Kunkle_hc.mg_go, aes(x = GO.name, y = score)) + 
  geom_point(aes(color = overlap), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](FDR))) +
  scale_y_continuous(limits = c(0, 9)) +
  theme(axis.text.y = element_text(size = 8)) + 
  scale_color_gradient(name = "genes", limits = c(0, 20), breaks = c(0, 5, 10, 15, 20))

ggplot(tau_hc.ec_go , aes(x = GO.name, y = score)) + 
  geom_point(aes(color = overlap), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](FDR))) +
  scale_y_continuous(limits = c(0, 7)) +
  theme(axis.text.y = element_text(size = 8)) +
  scale_color_gradient(name = "genes", limits = c(0, 20), breaks = c(0, 5, 10, 15, 20))

ggplot(tau_hc.ec2_go, aes(x = GO.name, y = score)) + 
  geom_point(aes(color = overlap), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](FDR))) +
  scale_y_continuous(limits = c(0, 2)) +
  theme(axis.text.y = element_text(size = 8)) + 
  scale_color_gradient(name = "genes", limits = c(0, 20), breaks = c(0, 5, 10, 15, 20))

#plot for KEGG enrichment
ggplot(tau_hc.ec_kegg, aes(x = KEGG.name, y = score)) + 
  geom_point(aes(color = overlap), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](FDR))) +
  scale_y_continuous(limits = c(0, 9)) +
  theme(axis.text.y = element_text(size = 8)) + 
  scale_color_gradient(name = "genes", limits = c(0, 20), breaks = c(0, 5, 10, 15, 20))

ggplot(tau_hc.ec2_kegg, aes(x = KEGG.name, y = score)) + 
  geom_point(aes(color = overlap), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](FDR))) +
  scale_y_continuous(limits = c(0, 9)) +
  theme(axis.text.y = element_text(size = 8)) + 
  scale_color_gradient(name = "genes", limits = c(0, 20), breaks = c(0, 5, 10, 15, 20))

ggplot(tau_hc.ca1_kegg, aes(x = KEGG.name, y = score)) + 
  geom_point(aes(color = overlap), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](FDR))) +
  scale_y_continuous(limits = c(0, 9)) +
  theme(axis.text.y = element_text(size = 8)) + 
  scale_color_gradient(name = "genes", limits = c(0, 20), breaks = c(0, 5, 10, 15, 20))
