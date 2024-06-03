#!/usr/bin/env Rscript
## This script plot individual-level heterozygosity estimates

## Load required packages
library(tidyverse)
library(cowplot)
## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
file <- args[2]
plot <- args[3]
sample_table_path <- args[4]
group_by <- args[5]
min_k <- args[6] %>% as.integer()
max_k <- args[7] %>% as.integer()

if(length(args)>7){
  pop_col <- args[8]
  pop <- args[9]
}

# indir <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/ohana/global"
# file <- "combined.subsetted"
# plot <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/figures/ohana/global/combined.subsetted.png"
# sample_table_path <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/docs/metadata.tsv"
# group_by <- "pca_cluster"
# min_k <- 2
# max_k <- 7 

# indir <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/ohana/local"
# file <- "combined.subsetted"
# plot <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/figures/ohana/local/vermilion.combined.subsetted.png"
# sample_table_path <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/docs/metadata.tsv"
# group_by <- "pca_cluster"
# min_k <- 2
# max_k <- 7 
# pop_col <- "common_name"
# pop <- "vermilion"

## read in sample table and subset by population if "local"
if(length(args)<=7){
  sample_table <- read_tsv(sample_table_path) %>%
    mutate(id_tmp = row_number() %>% as.character())
  indir_file <- str_c(indir, "/", file)
} else {
  sample_table <- read_tsv(sample_table_path) %>%
    filter(get(pop_col)==pop) %>%
    mutate(id_tmp = row_number() %>% as.character())
  indir_file <- str_c(indir, "/", pop, ".", file)
}
## read in admixture proportions
sample_size <- nrow(sample_table)
for (k in min_k:max_k) {
  genome_admix_k <-  read_tsv(paste0(indir_file, ".k", k,".q.matrix"), skip=1, col_names = FALSE) %>%
    as.matrix() %>%
    as_tibble() %>%
    bind_cols(sample_table) %>%
    pivot_longer(cols = 1:k, names_to = "anc", values_to = "p") %>%
    mutate(k = k)
  if (k == min_k) {
    genome_admix <- genome_admix_k
  } else {
    genome_admix <- bind_rows(genome_admix, genome_admix_k)
  }
}
## generate admixture plot
sample_order <- genome_admix %>%
  filter(k==max_k) %>%
  group_by(k, get(group_by), anc) %>%
  summarise(mean_p=mean(p)) %>%
  rename(!!group_by:=`get(group_by)`) %>%
  slice_max(order_by = mean_p) %>%
  ungroup() %>%
  semi_join(genome_admix, ., by=c(group_by, "anc", "k")) %>%
  arrange(get(group_by), -p) %>%
  .$id_tmp
admix_plot <- genome_admix %>%
  mutate(id_tmp=fct_relevel(id_tmp, sample_order)) %>%
  filter(p>0.001) %>%
  ggplot(aes(x = id_tmp, y = p, fill = anc, color = anc)) +
  geom_bar(stat = "identity", color="black", linewidth=0.3, width = 1) +
  #scale_fill_manual(values = color_palette) +
  #scale_color_manual(values = color_palette) +
  facet_grid((k) ~ get(group_by), scales = "free_x", space = "free_x", switch="y") +
  theme_cowplot() +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.line=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none",
        text = element_text(size=10),
        strip.text.x = element_text(size = 10, angle = 90),
        strip.text.y = element_text(size = 10, angle = 180))
ggsave(plot, admix_plot, width = sample_size*0.1, height = (max_k-min_k+2)*1, unit="cm")
