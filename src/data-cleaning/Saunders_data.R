#Preprocessing steps for Saunders et al data set

##### 1. load packages and data#######
### 1 load packages
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

#load functions for read dropviz data
source(here("src","tools","read_dropviz_data.R"))

##### 2 load data #######
all_files = list.files(here("raw","expr","Saunders"))
regions = c("Cerebellum_ALT","Cortex_noRep5_FRONTALonly","Cortex_noRep5_POSTERIORonly","EntoPeduncular","GlobusPallidus","Hippocampus","Striatum","SubstantiaNigra","Thalamus")
rg_formal = c("Cerebellum","Frontal Cortex","Posterior Cortex","Ento Peduncular","Globus Pallidus","Hippocampus","Striatum","Substantia Nigra","Thalamus")
tissues = c("CB","FC","PC","ENT","GP","HC","STR","SN","TH")
cell_type = c("Astrocytes","Endothelial","FibroblastLike","Microglia_Macrophage","Mural","Oligodendrocytes","Polydendrocytes")
##read in
dge = map(regions, ~loadSparseDge(here("raw","expr","Saunders",paste0("F_GRCm38.81.P60",.x,".raw.dge.txt.gz")))) %>% set_names(regions)
cell_info = map(regions, ~readRDS(here("raw","expr","Saunders",paste0("F_GRCm38.81.P60",.x,".cell_cluster_outcomes.RDS"))) %>% as_tibble(rownames="cellname")) %>% set_names(regions)
cell_type_voc = readRDS(here("raw","expr","Saunders","annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")) %>% as_tibble() 
non_neuron_ct = map(cell_type, ~readRDS(here("raw","expr","Saunders",all_files[grepl(all_files, pattern=.x)]))) %>% 
  set_names(cell_type) %>%
  map(~as_tibble(.x, rownames = "cellname")) %>% 
  map(~dplyr::rename(.x,non_neuron_cluster=subcluster)) 
non_neuron_ct = map2(non_neuron_ct, names(non_neuron_ct), ~mutate(.x, non_neuron_cluster = paste0(.y,"_",non_neuron_cluster))) %>%
  purrr::reduce(~rbind(.x,.y))

### 3 merge data
tot_genes = dge %>%  #keep only intersected genes
  map(~rownames(.x)) %>%
  purrr::reduce(~intersect(.x,.y))

data_expr = dge %>% 
  map(~.x[match(tot_genes, rownames(.x)),]) %>%
  map(~as(.x, "CsparseMatrix")) %>% 
  purrr::reduce(~cbind(.x,.y))

cell_info = cell_info %>%
  map2(tissues, ~mutate(.x, tissue_type = .y, subcluster=paste0(.y,"_",subcluster))) %>%
  purrr::reduce(~rbind(.x,.y)) 

# clean data annotation
# remove cells as doublet/outlier/small cells/without clear annotation
cell_info = cell_info %>% 
  mutate_if(is.factor, ~as.character(.)) %>%
  mutate(reason = ifelse(is.na(reason), "no_reason",reason)) %>%
  filter( reason=="curation"|reason == "no_reason") %>% #filter by annotation #1 # "small cell don't have annotation thus is also filtered out
  filter(!cellname %in% (non_neuron_ct %>% filter(!is.na(reason)) %>% pull(cellname))) %>%
  mutate(subcluster = ifelse( reason == "curation", paste0(subcluster,"-1"), subcluster))

# clean labels of neuron clusters
cell_type_voc = cell_type_voc %>%
  mutate(first_sub_marker = strsplit(type_marker, split="-") %>% map(~return(.[1])) %>% unlist ) %>%
  mutate(first_sub_marker= strsplit(first_sub_marker, split="[.-]") %>%  map(~return(.[1])) %>% unlist ) %>% 
  mutate(cluster=strsplit(subcluster, split="-") %>% map(~return(.[1])) %>% unlist) #split first subcluster by the primary marker

cluster_anno = rbind(data.frame(tissue = "FC", 
                                cluster = c(1,2,3,4,5,6,7,11),
                                cluster_anno = c("Cplx3+/Synpr+ interneurons", "MGE-derived cortical interneurons","Deep layer pyramidal cells","Deep layer pyramidal cells-Layer 5b","Claustrum","Superficial/Deep layer pyramidal cells","Deep layer pyramidal cells--layer 5","Cajal-Retzius")),
                     data.frame(tissue = "PC", 
                                cluster = c(1,2,3,4,5,6,7),
                                cluster_anno = c("Layer6_Subplate_Syt6", "Layer23_5_6_Subiculum","Layer5_Bcl6","CGE_Cplx3-Synpr interneurons","MGE_Sst-Pvalb interneurons","Layer5_Bcl6","RSG_Tshz2-Shisa8")),
                     data.frame(tissue = "HC", 
                                cluster = c(1,2,3,4,5,6,14),
                                cluster_anno = c("Gad2 interneurons", "Subiculum_Slc17a6","Subiculum_Entorhinal_Nxph3","Neuron_Dentate_C1ql2","Subiculum_Postsubiculum_Entorhinal","CA2_CA3","Cajal-Retzius")),
                     data.frame(tissue = "TH", 
                                cluster = c(1,2,3),
                                cluster_anno = c("Habenula_Tac2", "Rora","Gad2-Ahi1")),
                     data.frame(tissue = "STR", 
                                cluster = c(10,11,12,13,14,15),
                                cluster_anno = c("dSPN_Drd1", "iSPN_Drd1","cholinergic interneurons","SPN_BNST_Amygdala_Otof","Fast-spiking interneurons","Sst interneurons")),
                     data.frame(tissue = "GP", 
                                cluster = c(1,2,3),
                                cluster_anno = c("cholinergic interneurons", "LOT_GP_VGP_SI","SPN_BMA")),
                     data.frame(tissue = "ENT", 
                                cluster = 4,
                                cluster_anno = "Neuron_Syt1"),
                     data.frame(tissue = "SN", 
                                cluster = c(1,2,3,4),
                                cluster_anno = c("CA1_C1ql3","Thalamus","Gad2 neurons","DA neurons")),
                     data.frame(tissue="CB",
                                cluster = c(1,2,3,4),
                                cluster_anno = c("Granular neurons_Gabra6","Purkinje Neurons","Pvalb interneurons","Other interneurons"))) 

cell_type_voc = cell_type_voc %>% 
  left_join(cluster_anno %>% mutate(cluster=as.character(cluster)), by=c("tissue"="tissue","cluster"="cluster")) %>%
  mutate(cluster_anno = ifelse(is.na(cluster_anno) | class!="NEURON","",cluster_anno))

cell_type_voc = cell_type_voc %>%
  mutate(fine_cluster = ifelse(!tissue %in% c("CB","FC") , common_name, first_sub_marker))  %>% #interneurons are labelled by marker or 
  mutate(fine_cluster = recode(fine_cluster,
                              "Neurofilament, Layer4?" ="Neurofilament (possible Layer4)" ,
                              "Layer 2/3, IEG+" = "Layer 2/3",
                              "Layer 2/3, Cadm2+" = "Layer 2/3",
                              "Layer 5a, BC006965+" = "Layer 5a",
                              "Retrosplenial cortex (RSG), layer 5a" = "Layer 5a",
                              "Layer 5a, IEG+" = "Layer 5a",
                              "Cortical subplate interneuron" = "Cortical subplate interneuron, Cplx3",
                              "Neuron.Gad1Gad2-Slc17a8.Synpr-Sncg-Yjefn3" = "CGE_Synpr interneurons",
                              "Neuron.Gad1Gad2.Synpr-Sncg" = "CGE_Synpr interneurons",
                              "Neuron.Gad1Gad2.Synpr-Nnat" = "CGE_Synpr interneurons",
                              "Neuron.Gad1Gad2.Synpr-Crispld2" = "CGE_Synpr interneurons",
                              "Neuron.Gad1Gad2-Chat.Synpr-Slc5a7" = "CGE_Synpr interneurons",
                              "Neuron.Gad1Gad2.Synpr-Pcdh11x" = "CGE_Synpr interneurons",
                              "Neuron.Gad1Gad2.Pvalb-Unc5b" = "MGE_Pvalb interneurons",
                              "Neuron.Gad1Gad2.Pvalb-Nefm" = "MGE_Pvalb interneurons",
                              "Neuron.Gad1Gad2.Pvalb-Cartpt" = "MGE_Pvalb interneurons",
                              "Neuron.Gad1Gad2-Th.Nr4a2" = "MGE_Pvalb interneurons, Nr4a2",
                              "Neuron.Gad1Gad2.Sst-Crh" = "MGE_Sst interneurons",
                              "Neuron.Gad1Gad2.Sst-Pdyn" = "MGE_Sst interneurons",
                              "Neuron.Gad1Gad2.Sst-Nr2f2" = "MGE_Sst interneurons",
                              "Entorhinal cortex?" = "Possible entorhinal cortex",
                              "Layer 5, Retrosplenial cortex (RSG) and Subiculum, Tshz2+" = "RSG_Tshz2-Shisa8",
                              "Ventral Globus Pallidus Externus (vGP) / Ventral Pallidum (VP) / magnocellular preoptic nucleus (MCPO)"="ventral GP (vGP)",
                              "Substantia Innominata (SI), dorsal, bordering GP" = "Substantia Innominata (SI)",
                              "Striatum, iSPNs, Pde1c +"="Striatum, iSPNs",
                              "Interneuron, Neurogliaform1"="Interneuron, Neurogliaform", "Interneuron, Neurogliaform2"="Interneuron, Neurogliaform","Interneuron, Neurogliaform3"="Interneuron, Neurogliaform",
                              "Interneuron, OLM1"="Interneuron, OLM", 
                              "Interneuron, OLM2 (CA1 enriched?)"="Interneuron, OLM (possible CA1 enriched)", 
                              "Interneuron, OLM3 (CA1 enriched?)"="Interneuron, OLM (possible CA1 enriched)",
                              "Interneuron, Olm" = "Interneuron, OLM", 
                              "Interneuron, OLM4 (Dentate enriched?)"="Interneuron, OLM (possible Dentate enriched)",
                              "Entorhinal cortex (IEG)"="Entorhinal cortex",
                              "CA1 Principal cells (Anterior)"="CA1 Principal cells", 
                              "Medial entorrhinal cortexm, Reln+/Cbln1+"="Medial entorhinal cortex 1","Medial entorrhinal cortex, Cbln1+/Trps1+"="Medial entorhinal cortex 3",
                              "Medial entorrhinal cortex"="Medial entorhinal cortex 2",
                              "Lateral CA3 Principal cells"="CA3 Principal cells",
                              "eccentric SPN, iSPN-like markers, Th+" = "eccentric SPN, iSPN-like markers",
                              "iSPN, IEG+"="iSPN",
                              "dSPN, IEG+"="dSPN","dSPN, patch"="dSPN", "dSPN, lateral striatum"="dSPN",
                              "Fast-spiking interneuron, Pvalb+/Pnoc+"="Fast-spiking interneuron, Pvalb+", "Fast-spiking interneuron, Pvalb+/Rgs12+"="Fast-spiking interneuron, Pvalb+",
                              "Medial habenular, lateral portion"="Medial habenula",
                              " Supramammillary Nucleus (SuM) Dopaminergic?, Th+/Ddc+"="Supramammillary Nucleus (SuM)")) %>%#recode annotation and its typo
  mutate(fine_cluster = ifelse(grepl(pattern='CGE-derived ',x=fine_cluster), "Interneuron, candidate CGE-derived", fine_cluster)) %>% 
  mutate(fine_cluster = ifelse(tissue %in% c("CB","FC","PC") & grepl(pattern="[Ee]ntorhinal",x=common_name), paste0(fine_cluster," EC included"), fine_cluster)) %>%
  mutate(fine_cluster = ifelse(cluster_anno == "Granular neurons_Gabra6", "Granule cells",fine_cluster)) %>% 
  mutate(fine_cluster = ifelse(cluster_anno == "Purkinje Neurons", "Purkinje Neurons",fine_cluster)) %>% 
  mutate(fine_cluster = ifelse(tissue=="FC"&class=="NEURON",ifelse(grepl(pattern = "interneuron",x=cluster_anno),paste0(cluster_anno,".",fine_cluster), cluster_anno), fine_cluster)) %>%
  group_by(fine_cluster, class, tissue) %>% #if automatic cluster assignment failed 
  mutate(n=length(unique(cluster))) %>%
  mutate(fine_cluster=ifelse((class=="NEURON" & n>1), paste0(fine_cluster, ".cluster.",cluster),fine_cluster)) %>% #split neuron types (with same name but from differnt clusters)
  select_at(vars(-"n")) %>%
  ungroup %>%
  mutate(neuron_marker = strsplit(class_marker, split="-") %>% map(~.x[[1]]) %>% unlist) %>% mutate(subclass = ifelse(class!="NEURON",class, "")) %>%
  mutate(subclass = ifelse(subclass=="" & grepl("Gad1",x=neuron_marker), "GABAergic", subclass)) %>% 
  mutate(subclass = ifelse(subclass=="" & grepl("Slc17",x=neuron_marker), "Glutamatergic", subclass)) %>% 
  mutate(subclass = ifelse(subclass=="" & grepl("Chat",x=neuron_marker), "Cholinergic", subclass)) %>% 
  mutate(subclass = ifelse(subclass=="" & grepl("Th",x=neuron_marker), "Dopaminergic",subclass)) %>%
  group_by(fine_cluster,class,tissue) %>%
  mutate(n = length(unique(subclass))) %>% 
  mutate(fine_cluster = ifelse(class=="NEURON" & n>1, paste0(fine_cluster, ".", subclass), fine_cluster)) %>%
  select_at(vars(-"n")) %>%
  ungroup


#combine with non-neuron annotation
##rule:neurons are marked by "new cluster", all non neurons are marked by "cell cluster", but neurogenesis/mitotic + choroid plexus + ependyma cells are not markerd. We use original labels to markerd it
cell_info = cell_info %>%
  left_join( cell_type_voc %>% dplyr::select(class, tissue, tissue_subcluster, class_marker, type_marker, cluster_anno,common_name,full_name,subclass,fine_cluster), by=c("tissue_type"="tissue","subcluster"="tissue_subcluster")) %>% #add neuron annotation
  dplyr::select(cellname, tissue_type, class, cluster,subcluster,class_marker,type_marker, cluster_anno, full_name, common_name, subclass, fine_cluster ) %>%
  left_join(non_neuron_ct %>% distinct(cellname, non_neuron_cluster) %>% mutate(non_neuron_cluster = gsub(pattern="_1-[0-9]*",replacement = "", x=non_neuron_cluster)), by="cellname") %>% #merge non-neuron cell annotation
  mutate(fine_cluster = ifelse(class=="NEURON", paste0("Neurons_",fine_cluster), non_neuron_cluster)) %>% 
  mutate(fine_cluster=ifelse(class=="CHOROID_PLEXUS", "Choroid_plexus", fine_cluster)) %>% #annotate the cell type without non_neuron annotation
  mutate(fine_cluster=ifelse(class=="NEUROGENESIS" | class=="MITOTIC", "Neurogenesis_mitosis", fine_cluster)) %>%
  mutate(fine_cluster=ifelse(class=="EPENDYMAL", "Ependymal", fine_cluster)) %>%
  drop_na(fine_cluster) %>% #drop these non-neuronal cells without annotated class 
  mutate(class = ifelse(class=="NEURON","Neurons",fine_cluster)) %>% #modify class/subclass to make it the same as non_neuron annotation and fine clusters
  mutate(subclass = ifelse(class=="Neurons", subclass, class)) %>%
  mutate(cluster_anno = ifelse(cluster_anno=="", class, cluster_anno)) %>%  
  mutate(fine_cluster = paste0(tissue_type, ".",fine_cluster)) %>%
  dplyr::select(-non_neuron_cluster) %>%
  dplyr::rename(region = tissue_type)

#### 4 merge to sce
if (!require("SingleCellExperiment")){
  if (!requidplyr::renamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SingleCellExperiment")
  library("SingleCellExperiment")
} 
data_expr = data_expr[,match(cell_info$cellname, colnames(data_expr))]

brain_sce = SingleCellExperiment(list(counts=data_expr),
                                 colData = DataFrame(cell_info),
                                 rowData = DataFrame(symbol = rownames(data_expr)))


save(brain_sce, cell_type_voc, file=here("data","expr","Saunders","Saunders_clean.rda"))
