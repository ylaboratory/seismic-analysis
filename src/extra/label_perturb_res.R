library(tidyverse)
library(dplyr)
library(here)
library(SingleCellExperiment)
library(magrittr)

l1_similarity_normalized <- function(x, y, valid_idx = NULL) {
  if (!is.null(valid_idx)) {
    x <- x[valid_idx]
    y <- y[valid_idx]
  }
  valid_pos <- which(!is.na(x) & !is.na(y) & is.finite(x) & is.finite(y))
  x_valid <- x[valid_pos]
  y_valid <- y[valid_pos]
  
  
  x_norm <- (x_valid - mean(x_valid)) / sd(x_valid)
  y_norm <- (y_valid - mean(y_valid)) / sd(y_valid)
  
  l1_dist <- sum(abs(x_norm - y_norm))
  
  return(1-l1_dist/(sum(abs(x_norm)) + sum(abs(y_norm))))
}

#import the 
parameter_df <- read.table(here("data", "expr", "score_robustness", "cell_type_sample_by_tissue", "new_para_df.txt"), header = T, sep = "\t")

all_ref <- map(parameter_df$output_header, ~read.table(paste0(.x, ".all.original_score_df.txt"), header = T)) %>%
  map(~as_tibble(.x))

all_sample <- map(parameter_df$output_header, ~ {header = .x; map(1:5, ~read.table(paste0(header, ".all.sample_",.x,".score_df.txt"), header = T))})

#measuring the distance 
all_ref <- all_ref %>% 
  map(~ {data_tbl <- .x; map(colnames(data_tbl)[2:ncol(data_tbl)], ~dplyr::select(data_tbl, all_of(c("gene", .x))))}) %>%
  map(~set_names(.x, map(.x, ~colnames(.x)[2]) %>% unlist )) %>%
  map(~map(.x, ~set_colnames(.x, c("gene", "reference"))))

all_sample <- all_sample %>% 
  map(~map(.x, ~ {data_tbl <- .x; map(colnames(data_tbl)[2:ncol(data_tbl)], ~dplyr::select(data_tbl, all_of(c("gene", .x))))})) %>%
  map(~list_transpose(.x)) %>%
  map(~map(.x, ~purrr::reduce(.x, ~full_join(.x, .y, by = "gene")))) %>%
  map(~map(.x, ~set_colnames(.x, c("gene", paste0("sample_", 1:5)))))  %>%
  map(~set_names(.x, names(all_ref[[1]])))

all_results <- map2(all_ref, all_sample, ~map2(.x, .y, ~full_join(.x, .y, by = "gene")))

all_results <- all_results %>%
  map(~map(.x, function(data_tbl) {
    map(paste0("sample_", 1:5), function(sample_col) {
     l1_similarity_normalized(data_tbl$reference, data_tbl[[sample_col]])
    })
  }))

all_results <- map(all_results, ~map2(.x, names(.x) , 
                                      ~data.frame(sample = paste0("sample_", 1:5), value = unlist(.x)) %>% set_colnames(c("sample", .y)))) %>%
  list_transpose() %>%
  map(~map2(.x, parameter_df$cell_type_no, ~mutate(.x, cell_type = .y))) %>%
  map(~map2(.x, parameter_df$cell_ratio, ~mutate(.x, cell_ratio = .y))) %>%
  map(~purrr::reduce(.x, ~rbind(.x, .y))) %>%
  map2(names(.),  
       ~set_colnames(.x, c("sample", "l1_similarity", "cell_type", "cell_ratio")) %>% mutate(score_type = .y)) %>%
  purrr::reduce(~rbind(.x,.y))

all_results <- all_results %>% filter(score_type %in% c("seismic_score", "si", "spc_score", "de_score")) %>%
  mutate(score_type = factor(score_type, levels = c("seismic_score",  "de_score", "spc_score", "si")))

ggplot(all_results,
       aes(x = cell_ratio, y = l1_similarity, fill = score_type)) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.alpha = 0.5) +
  theme_classic() +
  ggsci::scale_fill_npg(labels = c("seismic specificity", "DE score","Bryois gene specificity", "specificity index")) +
  xlab("ratio of extra sampled cells (relative to the original number of cells)") +
  ylab("normalized L1 similarity") +
  scale_x_discrete(labels=c("cr_0.05" = "5%", "cr_0.1" = "10%", "cr_0.25" = "25%", "cr_0.5" = "50%", "cr_0.75" = "75%", "cr_1" = "100%"))


###correlation across cell ratio
my_cor_func <- function(x, y) {
  valid_idx <- which(!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y))
  return(cor(x[valid_idx], y[valid_idx], method = "pearson"))
}

##measure stability
all_sample_t <- list_transpose(all_sample) 

cor_across_sample <- all_sample_t %>%
  map(~map(.x, ~ proxy::simil(t(.x[,2:6]), method = my_cor_func))) %>%
  map(~map(.x, ~.x %>% as.matrix() %>% as.data.frame() %>% as_tibble(rownames = "sample") %>% pivot_longer(cols = !sample, names_to = "sample2", values_to = "correlation") )) %>%
  map(~map(.x, ~mutate(.x, sample = factor(sample, levels = paste0("sample_",1:5)), sample2 = factor(sample2, levels = paste0("sample_",1:5))))) %>%
  map(~map(.x, ~filter(.x, as.numeric(sample) < as.numeric(sample2)))) %>%
  map(~map2(.x, parameter_df$cell_type_no, ~mutate(.x, cell_type = .y))) %>%
  map(~map2(.x, parameter_df$cell_ratio, ~mutate(.x, cell_ratio = .y))) %>%
  map(~purrr::reduce(.x, ~rbind(.x, .y))) %>%
  map2(names(.), ~mutate(.x, score_type = .y)) %>%
  purrr::reduce(~rbind(.x, .y))
  
ggplot(cor_across_sample %>% filter(score_type %in% c("seismic_score", "si", "spc_score", "de_score")),
       aes(x = cell_ratio, y = correlation, fill = score_type)) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.alpha = 0.5)
