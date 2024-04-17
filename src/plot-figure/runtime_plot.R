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
color_mapping_vec = c("scDRS" ="#46B8DAFF" , "seismic" = "#D43F3AFF" , "S-MAGMA"="#EEA236FF" , "FUMA"="#5CB85CFF")

#load data
all_time = list("seismic" = here("results","runtime","seismic","runtime.txt"),
                "scDRS" = here("results","runtime","scDRS","runtime.txt"),
                "FUMA" = here("results","runtime","FUMA","runtime.txt"),
                "S-MAGMA" = here("results","runtime","S-MAGMA","runtime.txt")) %>%
  map(~read.table(.x,sep="\t",header=T)) %>% 
  map(~as_tibble(.x)) %>% 
  map(~group_by(.x, size)) %>% 
  map(~mutate(.x, total_time = processing_time + association_time)) %>%
  map(~mutate(.x, time_sd = sd(total_time),total_time = mean(total_time))) %>% 
  map(~ungroup(.x)) %>% 
  map(~distinct(.x, size,total_time,time_sd)) %>%
  map2(names(.), ~mutate(.x, method = .y)) %>%
  purrr::reduce(~rbind(.x, .y)) %>%
  mutate(size = as.numeric(gsub(pattern="k",replacement = "", x=size)))

#total plot
ggplot(all_time, aes(x = size, y = total_time, group = method, color = method)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = total_time - time_sd, ymax = total_time + time_sd), width = 0.02, alpha=0.7) +
  scale_y_log10()+ 
  scale_x_log10() + 
  labs(title = "Total Time by Data Size and Method",
       x = "number of cells (thousands, log scale)",
       y = "total running time (seconds, log scale)") +
  scale_color_manual(values=color_mapping_vec)+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#print
all_time %>% write.csv(here("results","runtime","runtime_all.csv"), row.names = F)
