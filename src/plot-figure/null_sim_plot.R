#plot for the simulation results
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

#quantiles 
qqpoint <- floor(10000/10^ seq(0,3,0.1))
qqpos <- -log10(qqpoint/10000)


##### plot for scrambled expression #####
#load results
# null_res_expr = list.files(here("results","null_sim","expr_rs"),pattern = "\\.null_res\\.txt") %>%
#  set_names(gsub(pattern = "\\.null_res\\.txt", replacement = "", x=unlist(.))) %>%
#  map(~read.table(here("results","null_sim","expr_rs",.x),header = T)) %>% 
null_res_expr <- list.files(here("results","null_sim","expr_rs"),pattern = "new_ds_[0-9]*.null_res\\.txt",full.names = T) %>%
  set_names(str_extract(string = ., pattern = "ds_[0-9]*")) %>%
  map(~read.table(.x,header = T)) %>% 
  map2(names(.), ~mutate(.x, ds = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% 
  set_colnames(c("index","seismic","S-MAGMA","FUMA","ds"))

null_res_expr_scdrs <- list.files(here("results","null_sim","expr_rs"),pattern = "\\.scdrs_null_res\\.txt") %>%
  set_names(gsub(pattern = "\\.scdrs_null_res\\.txt", replacement = "", x=unlist(.))) %>%
  map(~read.table(here("results","null_sim","expr_rs",.x),header = T)) %>% 
  map(~mutate(.x, scDRS = pnorm(assoc_mcz, lower.tail = F)))  %>%
  map(~select(.x, -assoc_mcz)) %>%
  map2(names(.), ~mutate(.x, ds = .y)) %>%
  purrr::reduce(~rbind(.x,.y))

null_res_expr_all <- null_res_expr %>% 
  left_join(null_res_expr_scdrs, by=c("index"="index", "ds" = "ds")) %>% 
  select(any_of(c("index", "ds", "seismic", "scDRS", "S-MAGMA", "FUMA"))) 

null_res_expr_sum <- map(c("seismic","scDRS", "S-MAGMA", "FUMA"), ~select(null_res_expr_all, all_of(c("index","ds",.x)))) %>% #summary info
  set_names(c("seismic","scDRS", "S-MAGMA", "FUMA")) %>%
  map2(names(.), ~pivot_wider(.x, names_from = "ds", values_from = .y)) %>% 
  map(~mutate(.x, across( matches("ds"), ~base::sort(.,na.last=T)))) %>%
  map(~dplyr::slice(.x, qqpoint)) %>% #slice quantiles
  map(~rowwise(.x)) %>% 
  map(~mutate_if(.x,is.double, ~-log10(.))) %>% #-log10
  map(~mutate(.x,mean =  mean(c_across(starts_with("ds_"))), sd = sd(c_across(starts_with("ds_"))) )) %>%  #sd and mean
  map(~ungroup(.x)) %>%
  map(~mutate(.x,exp_q = qqpos)) %>%
  map(~mutate(.x,mean_max = mean+sd)) %>% 
  map(~mutate(.x,mean_min = mean-sd)) %>%
  map(~select(.x, mean, exp_q, mean_max,mean_min)) %>%
  map2(names(.), ~mutate(.x, method = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% #combine together
  rowwise() %>%
  mutate(mean_min = max(0,mean_min)) %>%
  mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))) 
  
ggplot(data = null_res_expr_sum, aes(x=exp_q, y=mean, ymin = mean_min, ymax=mean_max, fill=method)) +  
  geom_line(aes(colour = method)) +
  geom_abline(slope=1,linetype="dashed") + 
  geom_ribbon(alpha=0.5) + 
  xlab("expected -log10(P-value) quantiles") + 
  ylab("observed -log10(P-value) quantiles") + 
  facet_wrap(vars(method),nrow=2) +
  scale_fill_manual(values=color_mapping_vec)+
  scale_color_manual(values=color_mapping_vec)+
  ggtitle("Null distribution of P-value between traits and cell types, shuffled gene index") +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5,size=12),axis.text.y = element_text(size=8), strip.background = element_blank())

null_res_expr_sum %>% write.csv(here("results","null_sim","null_expr_all.csv"), row.names = F)

##### plot for scrambled gene set #####
null_res_gs <- list.files(here("results","null_sim","gs_rs"),pattern = "new_ds_[0-9]*.null_res\\.txt",full.names = T) %>%
  set_names(str_extract(string = ., pattern = "ds_[0-9]*")) %>%
  map(~read.table(.x,header = T)) %>%
  map2(names(.), ~mutate(.x, ds = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% 
  set_colnames(c("index","seismic","ds"))

null_res_gs_scdrs <- list.files(here("results","null_sim","gs_rs"),pattern = "\\.scdrs_null_res\\.txt") %>%
  set_names(gsub(pattern = "\\.scdrs_null_res\\.txt", replacement = "", x=unlist(.))) %>%
  map(~read.table(here("results","null_sim","gs_rs",.x),header = T)) %>% 
  map(~mutate(.x, scDRS = pnorm(assoc_mcz, lower.tail = F)))  %>%
  map(~select(.x, -assoc_mcz)) %>%
  map2(names(.), ~mutate(.x, ds = .y)) %>%
  purrr::reduce(~rbind(.x,.y))

null_res_gs_all <- null_res_gs %>% 
  left_join(null_res_gs_scdrs, by=c("index"="index", "ds" = "ds")) %>% 
  select(any_of(c("index", "ds", "seismic", "scDRS"))) 

null_res_gs_sum <- map(c("seismic","scDRS"), ~select(null_res_gs_all, all_of(c("index","ds",.x)))) %>% #summary info
  set_names(c("seismic","scDRS")) %>%
  map2(names(.), ~pivot_wider(.x, names_from = "ds", values_from = .y)) %>% 
  map(~mutate(.x, across( matches("ds"), ~base::sort(.,na.last=T)))) %>%
  map(~dplyr::slice(.x, qqpoint)) %>% #slice quantiles
  map(~rowwise(.x)) %>% 
  map(~mutate_if(.x,is.double, ~-log10(.))) %>% #-log10
  map(~mutate(.x,mean =  mean(c_across(starts_with("ds_"))), sd = sd(c_across(starts_with("ds_"))) )) %>%  #sd and mean
  map(~ungroup(.x)) %>%
  map(~mutate(.x,exp_q = qqpos)) %>%
  map(~mutate(.x,mean_max = mean+sd)) %>% 
  map(~mutate(.x,mean_min = mean-sd)) %>%
  map(~select(.x, mean, exp_q, mean_max,mean_min)) %>%
  map2(names(.), ~mutate(.x, method = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% #combine together
  rowwise() %>%
  mutate(mean_min = max(0,mean_min)) %>%
  mutate(method = factor(method, levels = c("seismic","scDRS"))) 

ggplot(data = null_res_gs_sum, aes(x=exp_q, y=mean, ymin = mean_min, ymax=mean_max, fill=method)) +  
  geom_line(aes(colour = method)) +
  geom_abline(slope=1,linetype="dashed") + 
  geom_ribbon(alpha=0.5) + 
  xlab("expected -log10(P-value) quantiles") + 
  ylab("observed -log10(P-value) quantiles") + 
  facet_wrap(vars(method),nrow=1) +
  scale_fill_manual(values=color_mapping_vec)+
  scale_color_manual(values=color_mapping_vec)+
  ggtitle("Null distribution of P-value between traits and cell types, shuffled gene index") +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5,size=12),axis.text.y = element_text(size=8), strip.background = element_blank())

null_res_gs_sum %>% write.csv(here("results","null_sim","null_gs_all.csv"), row.names = F)


##### plot for scrambled expression #####
#load results
# null_res_expr = list.files(here("results","null_sim","expr_rs"),pattern = "\\.null_res\\.txt") %>%
#  set_names(gsub(pattern = "\\.null_res\\.txt", replacement = "", x=unlist(.))) %>%
#  map(~read.table(here("results","null_sim","expr_rs",.x),header = T)) %>% 
null_res_expr <- list.files(here("results","null_sim","expr_rs"),pattern = "new_ds_[0-9]*.null_res\\.txt",full.names = T) %>%
  set_names(str_extract(string = ., pattern = "ds_[0-9]*")) %>%
  map(~read.table(.x,header = T)) %>% 
  map2(names(.), ~mutate(.x, ds = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% 
  set_colnames(c("index","seismic","S-MAGMA","FUMA","ds"))

null_res_expr_scdrs <- list.files(here("results","null_sim","expr_rs"),pattern = "\\.scdrs_null_res\\.txt") %>%
  set_names(gsub(pattern = "\\.scdrs_null_res\\.txt", replacement = "", x=unlist(.))) %>%
  map(~read.table(here("results","null_sim","expr_rs",.x),header = T)) %>% 
  map(~mutate(.x, scDRS = pnorm(assoc_mcz, lower.tail = F)))  %>%
  map(~select(.x, -assoc_mcz)) %>%
  map2(names(.), ~mutate(.x, ds = .y)) %>%
  purrr::reduce(~rbind(.x,.y))

null_res_expr_all <- null_res_expr %>% 
  left_join(null_res_expr_scdrs, by=c("index"="index", "ds" = "ds")) %>% 
  select(any_of(c("index", "ds", "seismic", "scDRS", "S-MAGMA", "FUMA"))) 

null_res_expr_sum <- map(c("seismic","scDRS", "S-MAGMA", "FUMA"), ~select(null_res_expr_all, all_of(c("index","ds",.x)))) %>% #summary info
  set_names(c("seismic","scDRS", "S-MAGMA", "FUMA")) %>%
  map2(names(.), ~pivot_wider(.x, names_from = "ds", values_from = .y)) %>% 
  map(~mutate(.x, across( matches("ds"), ~base::sort(.,na.last=T)))) %>%
  map(~dplyr::slice(.x, qqpoint)) %>% #slice quantiles
  map(~rowwise(.x)) %>% 
  map(~mutate_if(.x,is.double, ~-log10(.))) %>% #-log10
  map(~mutate(.x,mean =  mean(c_across(starts_with("ds_"))), sd = sd(c_across(starts_with("ds_"))) )) %>%  #sd and mean
  map(~ungroup(.x)) %>%
  map(~mutate(.x,exp_q = qqpos)) %>%
  map(~mutate(.x,mean_max = mean+sd)) %>% 
  map(~mutate(.x,mean_min = mean-sd)) %>%
  map(~select(.x, mean, exp_q, mean_max,mean_min)) %>%
  map2(names(.), ~mutate(.x, method = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% #combine together
  rowwise() %>%
  mutate(mean_min = max(0,mean_min)) %>%
  mutate(method = factor(method, levels = c("seismic","scDRS","FUMA","S-MAGMA"))) 

ggplot(data = null_res_expr_sum, aes(x=exp_q, y=mean, ymin = mean_min, ymax=mean_max, fill=method)) +  
  geom_line(aes(colour = method)) +
  geom_abline(slope=1,linetype="dashed") + 
  geom_ribbon(alpha=0.5) + 
  xlab("expected -log10(P-value) quantiles") + 
  ylab("observed -log10(P-value) quantiles") + 
  facet_wrap(vars(method),nrow=2) +
  scale_fill_manual(values=color_mapping_vec)+
  scale_color_manual(values=color_mapping_vec)+
  ggtitle("Null distribution of P-value between traits and cell types, shuffled gene index") +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5,size=12),axis.text.y = element_text(size=8), strip.background = element_blank())

null_res_expr_sum %>% write.csv(here("results","null_sim","null_expr_all.csv"), row.names = F)

##### plot for new fuma #####
fuma_null_res_gs <- list.files(here("results","null_sim","expr_rs"),pattern = "fuma_ds_[0-9]*.null_res\\.txt",full.names = T) %>%
  set_names(str_extract(string = ., pattern = "ds_[0-9]*")) %>%
  map(~read.table(.x,header = T)) %>%
  map2(names(.), ~mutate(.x, ds = .y)) %>%
  purrr::reduce(~rbind(.x,.y)) %>% 
  set_colnames(c("index","FUMA","ds"))

fuma_null_res_sum <- fuma_null_res_gs %>%
  pivot_wider( names_from = "ds", values_from = FUMA) %>% 
  mutate( across( matches("ds"), ~base::sort(.,na.last=T))) %>%
  dplyr::slice( qqpoint) %>% #slice quantiles
  rowwise() %>% 
  mutate_if(is.double, ~-log10(.)) %>% #-log10
  mutate(mean =  mean(c_across(starts_with("ds_"))), sd = sd(c_across(starts_with("ds_"))) ) %>%  #sd and mean
  ungroup() %>%
  mutate(exp_q = qqpos) %>%
  mutate(mean_max = mean+sd) %>% 
  mutate(mean_min = mean-sd) %>%
  select( mean, exp_q, mean_max,mean_min) %>%
  rowwise() %>%
  mutate(mean_min = max(0,mean_min)) 

ggplot(data = fuma_null_res_sum, aes(x=exp_q, y=mean, ymin = mean_min, ymax=mean_max)) +  
  geom_line(color=color_mapping_vec["FUMA"]) +
  geom_abline(slope=1,linetype="dashed") + 
  geom_ribbon(alpha=0.5,fill=color_mapping_vec["FUMA"]) + 
  xlab("expected -log10(P-value) quantiles") + 
  ylab("observed -log10(P-value) quantiles") + 
  xlim(0,3) + ylim(0,3) +
  ggtitle("Null distribution of P-value between traits and cell types, shuffled expression") +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5,size=12),axis.text.y = element_text(size=8), strip.background = element_blank())
