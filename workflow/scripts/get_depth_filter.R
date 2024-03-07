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

## read in empirical depth distribution
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
## find the mode of the distribution
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
# subset the chunk of distribution between 50% and 150% of the mode
depth_lower_bound <- depth_mode*0.5
depth_upper_bound <- depth_mode*1.5
depth_hist_subset <- depth_hist %>%
  filter(depth > depth_lower_bound, depth < depth_upper_bound)
n_sites_subset <- depth_hist_subset %>%
  pull(count) %>%
  sum()
## fit a truncated normal distribution to the subsetted chunk
set.seed(42)
depth_dist_subset <- depth_hist_subset %>%
  uncount(count) %>%
  pull(depth) %>%
  sample(100000)
fitted_dist <- depth_dist_subset %>%
  fitdist("tnorm", start = list(mean=depth_mode, sd=depth_mode/5), fix.arg = list(a = depth_lower_bound, b= depth_upper_bound), discrete = TRUE)
## use numbers of sd away from the mean as depth filters
fitted_mean <- fitted_dist$estimate[1]
fitted_sd <- fitted_dist$estimate[2]
min_filter <- round(fitted_mean-fitted_sd*n_sd)
max_filter <- round(fitted_mean+fitted_sd*n_sd)
fitted_hist <- tibble(depth=(depth_lower_bound+1):(depth_upper_bound-1)) %>%
  mutate(count=dtnorm(depth, mean=fitted_mean, sd = fitted_sd, a = depth_lower_bound, b= depth_upper_bound)) %>%
  mutate(count=count*n_sites_subset)
## plot depth distribution and depth filters
depth_plot <- depth_hist %>%
  filter(count <= count_mode * 1.5) %>%
  ggplot(aes(x=depth, y=count)) +
  geom_line(data=fitted_hist, color="blue", linewidth=1) +
  geom_line() +
  geom_vline(xintercept = c(min_filter, max_filter), color = "red") +
  annotate("text", x= c(min_filter, max_filter), y = 0, label=c(min_filter, max_filter), hjust = -0.2, vjust = 0, color="red") +
  annotate("text", Inf, Inf, label=str_c("fitted_mean=", round(fitted_mean,1)), hjust = 1.1, vjust= 1.1, color="blue") +
  annotate("text", Inf, Inf, label=str_c("fitted_sd=", round(fitted_sd,1)), hjust = 1.1, vjust= 3.1, color="blue") +
  xlim(c(NA, fitted_mean*2)) +
  cowplot::theme_minimal_grid()
## output chosen depth filters as a text file
tibble(min_filter, max_filter, fitted_mean, fitted_sd) %>%
  write_tsv(str_c(indir, "/depth_filter.tsv"))
## output the plot
ggsave(str_c(plot_dir, "/depth_filter.png"), depth_plot, width = 6, height = 3, units = "in")
