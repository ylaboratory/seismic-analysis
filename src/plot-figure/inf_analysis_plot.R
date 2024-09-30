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
  as_tibble()
Kunkle_hc.mg_go <- read.table(here("results","Saunders","inf_analysis","Kunkle_hc_micro_go.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(p.adjust)) %>%
  head(n=10)

tau_hc.ec_dfbetas <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec_dfbetas.txt"), header = T, sep = "\t") %>%
  as_tibble()
tau_hc.ec_go <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec_go.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(p.adjust)) %>% 
  head(n=10)

tau_hc.ec2_dfbetas <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec2_dfbetas.txt"), header = T, sep = "\t") %>%
  as_tibble()
tau_hc.ec2_go <- read.table(here("results","Saunders","inf_analysis","tau_hc_ec2_go.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(score = -log10(p.adjust)) %>% 
  head(n=10)

#plot for influential genes
seismicGWAS::plot_inf_genes(Kunkle_hc.micro_dfbetas , gene_col = "genename")

seismicGWAS::plot_inf_genes(tau_hc.ec2_dfbetas, gene_col = "genename")

seismicGWAS::plot_inf_genes(tau_hc.ec_dfbetas, gene_col = "genename")

#plot for GO enrichment
ggplot(Kunkle_hc.mg_go, aes(x = Description, y = score)) + 
  geom_point(aes(color = Count), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](Pvalue))) +
  scale_y_continuous(limits = c(0, 9)) +
  theme(axis.text.y = element_text(size = 8)) + 
  guides(color=guide_legend(title="gene"))

ggplot(tau_hc.ec_go , aes(x = Description, y = score)) + 
  geom_point(aes(color = Count), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](Pvalue))) +
  scale_y_continuous(limits = c(0, 7)) +
  theme(axis.text.y = element_text(size = 8)) + 
  guides(color=guide_legend(title="gene"))

ggplot(tau_hc.ec2_go, aes(x = Description, y = score)) + 
  geom_point(aes(color = Count), size = 5) +
  coord_flip() + 
  theme_classic() +
  labs(x = "", y = expression(-log[10](Pvalue))) +
  scale_y_continuous(limits = c(0, 2)) +
  theme(axis.text.y = element_text(size = 8)) + 
  guides(color=guide_legend(title="gene"))

