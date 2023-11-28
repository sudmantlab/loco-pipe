#!/usr/bin/env Rscript
## This script plot individual-level heterozygosity estimates

## Load required packages
library(tidyverse)
library(cowplot)
## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
cov <- args[1]
plot <- args[2]
sample_table_path <- args[3]
color_by <- args[4]
if(length(args)>4){
  pop_col <- args[5]
  pop <- args[6]
}

# cov <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/pcangsd/local/vermilion.combined.subsetted.cov"
# plot <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/figures/pcangsd/local/vermilion.combined.subsetted.png"
# sample_table_path <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/docs/metadata.tsv"
# color_by <- "pca_cluster"
# pop_col <- "common_name"
# pop <- "vermilion"

if(length(args)<=4){
  sample_table <- read_tsv(sample_table_path)
} else {
  sample_table <- read_tsv(sample_table_path) %>%
    filter(get(pop_col)==pop)
}
sample_size <- nrow(sample_table)
c <- read_delim(cov, col_names = FALSE, delim=" ") %>%
  as.matrix()
e <- eigen(c)
e_values <- e$values
var_explained <- round(e_values/sum(e_values)*100, 1)
fix_names <- function(x) str_c("PC", seq_along(x))
e_vectors <- as_tibble(e$vectors, .name_repair = fix_names) 
pca_table <- sample_table %>%
  bind_cols(e_vectors) 
## note: we should subset less stringently
axis_1 <- 1
axis_2 <- 2
plot_pca <-function(axis_1, axis_2){
  pca_table %>%
    ggplot(aes(x=get(str_c("PC", axis_1)), 
               y=get(str_c("PC", axis_2)),
               color=get(color_by))) +
    geom_point(shape=21) +
    scale_color_discrete(name = color_by) +
    labs(x=str_c("PC", axis_1, " (", var_explained[axis_1],"%)"),
         y=str_c("PC", axis_2, " (", var_explained[axis_2],"%)")) +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.ticks = element_blank(),
          axis.text = element_blank())
}
legend <- get_legend(plot_pca(1,2) + theme(legend.position = "right"))
if (sample_size < 4) {
  (plot_pca(1,2) + theme(legend.position = "right")) %>% 
    ggsave(plot=., filename = plot, width = 4, height = 2, units = "in")
} else if (sample_size < 6) {
  plot_grid(plot_grid(plot_pca(1,2), plot_pca(3,4), ncol = 1),
            legend,
            nrow = 1,
            rel_widths = c(1, 1)
  ) %>% 
    ggsave(plot=., filename = plot, width = 5, height = 4, units = "in")
} else if (sample_size < 8) {
  plot_grid(plot_grid(plot_pca(1,2), plot_pca(3,4), ncol = 1),
            plot_grid(plot_pca(5,6), NULL, ncol = 1),
            legend,
            nrow = 1,
            rel_widths = c(1, 1, 1)
  ) %>% 
    ggsave(plot=., filename = plot, width = 7.5, height = 4, units = "in")
} else {
  plot_grid(plot_grid(plot_pca(1,2), plot_pca(3,4), ncol = 1),
            plot_grid(plot_pca(5,6), plot_pca(7,8), ncol = 1),
            legend,
            nrow = 1,
            rel_widths = c(1, 1, 1)
  ) %>% 
    ggsave(plot=., filename = plot, width = 7.5, height = 4, units = "in")
}
