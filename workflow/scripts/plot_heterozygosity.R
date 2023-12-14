#!/usr/bin/env Rscript
## This script plot individual-level heterozygosity estimates

## Load required packages
library(tidyverse)
library(cowplot)
## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
outdir <- args[2]
sample_table_path <- args[3]
color_by <- args[4]

# indir <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/angsd/heterozygosity"
# outdir <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/figures/heterozygosity"
# sample_table_path <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/docs/metadata.tsv"
# color_by <- "pca_cluster"

## read in the data
sample_table <- read_tsv(sample_table_path)
het <- NULL
for (i in 1:nrow(sample_table)){
  sfs_i <- read_delim(str_c(indir, "/", sample_table$sample_name[i], ".sfs"), col_names = FALSE) %>%
    as.matrix() %>%
    as.vector()
  het_i <- sfs_i[2]/sum(sfs_i, na.rm = TRUE)
  het <- c(het, het_i)
}

## plot a boxplot
het_plot <- mutate(sample_table, het=het) %>%
  ggplot(aes(x=fct_rev(get(color_by)), y=het, fill=fct_rev(get(color_by)))) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0) +
  labs(x=color_by, y="heterozygosity") +
  coord_flip() +
  theme_cowplot() +
  theme(legend.position = "none")

n_row <- sample_table %>% pull({color_by}) %>% unique() %>% length()
ggsave(str_c(outdir, "/heterozygosity.png"), het_plot, width = 6, height = n_row*1, units = "in")
