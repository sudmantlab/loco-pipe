#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
plot_dir <- args[2]
chr_table <- args[3]
n_sd <- as.double(args[4])

library(tidyverse)
library(cowplot)
library(fitdistrplus)
library(extraDistr)
#indir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/angsd/get_depth_global/"
#chr_list <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/docs/chr_list.txt"
chrs <- read_tsv(chr_table, col_names = FALSE) %>% pull(1)
for (i in seq_along(chrs)){
  depth_count_tmp <- read_lines(str_c(indir, "/", chrs[i], ".depthGlobal")) %>%
    str_split(pattern = "\t") %>%
    .[[1]] %>%
    .[1:10001] %>%
    as.integer()
  if (i == 1){
    depth_count <- depth_count_tmp
  } else {
    depth_count <- depth_count + depth_count_tmp
  }
}
depth_hist <- depth_count %>%
  tibble(depth=0:10000, count=.) %>%
  filter(depth > 0)
depth_quantile <- depth_hist %>%
  uncount(count) %>%
  pull(depth) %>%
  quantile(probs = 0.25)
depth_mode <- depth_hist %>%
  filter(depth>depth_quantile) %>%
  slice_max(count) %>%
  pull(depth)
count_mode <- depth_hist %>%
  filter(depth>depth_quantile) %>%
  slice_max(count) %>%
  pull(count)
depth_lower_bound <- depth_mode*0.5
depth_upper_bound <- depth_mode*1.5
depth_hist_subset <- depth_hist %>%
  filter(depth > depth_lower_bound, depth < depth_upper_bound)
n_sites_subset <- depth_hist_subset %>%
  pull(count) %>%
  sum()
set.seed(42)
depth_dist_subset <- depth_hist_subset %>%
  uncount(count) %>%
  pull(depth) %>%
  sample(100000)
fitted_dist <- depth_dist_subset %>%
  fitdist("tnorm", start = list(mean=depth_mode, sd=depth_mode/5), fix.arg = list(a = depth_lower_bound, b= depth_upper_bound), discrete = TRUE)
fitted_mean <- fitted_dist$estimate[1]
fitted_sd <- fitted_dist$estimate[2]
min_filter <- round(fitted_mean-fitted_sd*n_sd)
max_filter <- round(fitted_mean+fitted_sd*n_sd)
fitted_hist <- tibble(depth=(depth_lower_bound+1):(depth_upper_bound-1)) %>%
  mutate(count=dtnorm(depth, mean=fitted_mean, sd = fitted_sd, a = depth_lower_bound, b= depth_upper_bound)) %>%
  #dnbinom(x = 1:10000, size=fit$estimate[1], mu = fit$estimate[2]) %>% 
  mutate(count=count*n_sites_subset)
depth_plot <- depth_hist %>%
  filter(count <= count_mode * 1.5) %>%
  ggplot(aes(x=depth, y=count)) +
  geom_line(data=fitted_hist, color="blue", linewidth=1) +
  geom_line() +
  geom_vline(xintercept = c(min_filter, max_filter), color = "red") +
  annotate("text", x= c(min_filter, max_filter), y = Inf, label=c(min_filter, max_filter), hjust = -0.2, vjust = 3, color="red") +
  annotate("text", Inf, -Inf, label=str_c("fitted_mean=", round(fitted_mean,1)), hjust = 1.1, vjust= -20, color="blue") +
  annotate("text", Inf, -Inf, label=str_c("fitted_sd=", round(fitted_sd,1)), hjust = 1.1, vjust= -18, color="blue") +
  xlim(c(NA, fitted_mean*2)) +
  cowplot::theme_minimal_grid()
tibble(min_filter, max_filter, fitted_mean, fitted_sd) %>%
  write_tsv(str_c(indir, "/depth_filter.tsv"))
ggsave(str_c(plot_dir, "/depth_filter.png"), depth_plot)
