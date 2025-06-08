#plot for the simulation results
##### load packages and results ####
if (!require("here")) {
  install.packages("here")
  library("here")
}
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}


#color
color_mapping_vec = c("scDRS" ="#46B8DAFF" , "seismic" = "#D43F3AFF" , "S-MAGMA"="#EEA236FF" , "FUMA"="#5CB85CFF")

##### load runtime data #####
#function
load_and_extract <- function(data_set) {
  load(data_set)
  data.frame(processing_time = unlist(processing_time), association_time = unlist(group_time))
}

#parameters
sample_size <- c("10k", "25k", "50k", "100k", "150k", "200k", "250k", "300k", "400k", "500k")
ds <- paste0("ds_",1:5)

#load data
all_time_seismic <- expand.grid(size = sample_size, ds_index = ds) %>%
  rowwise() %>%
  mutate(data = map2(size, ds_index, ~load_and_extract(here("results", "runtime", .x, .y, paste0(.y, ".seismic.rda"))))) %>%
  unnest(data) %>%
  mutate(total_time = processing_time + association_time)

all_time_fuma <- expand.grid(size = sample_size, ds_index = ds) %>%
  rowwise() %>%
  mutate(data = map2(size, ds_index, ~load_and_extract(here("results", "runtime", .x, .y, paste0(.y, ".fuma.rda"))))) %>%
  unnest(data) %>% 
  mutate(total_time = processing_time + association_time)

all_time_magma <- expand.grid(size = sample_size, ds_index = ds) %>%
  rowwise() %>%
  mutate(data = map2(size, ds_index, ~load_and_extract(here("results", "runtime", .x, .y, paste0(.y, ".magma.rda"))))) %>%
  unnest(data) %>% 
  mutate(total_time = processing_time + association_time)

all_time_scdrs <- read.table(here("results","runtime","scDRS","runtime.txt"), header = T, sep="\t") %>%
  as_tibble() %>%
  mutate(total_time = processing_time+score_time+association_time)

#load data
all_time <- list("seismic" = all_time_seismic,
                 "scDRS" = all_time_scdrs,
                 "FUMA" = all_time_fuma,
                 "S-MAGMA" = all_time_magma) %>%
  map(~group_by(.x ,size)) %>%
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
  scale_y_continuous(limits = c(0,12000)) +
  scale_y_log10()+ 
  labs(title = "Total Time by Data Size and Method",
       x = "number of cells",
       y = "total running time (seconds, log scale)") +
  scale_color_manual(values=color_mapping_vec)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

#print
all_time %>% write.csv(here("results","runtime","runtime_all.csv"), row.names = F)
