#!/usr/bin/env Rscript
## This script read in PCA result in all windows and runs lostruct

## Load required packages
library(tidyverse)

## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
chr_table_path <- args[1]
in_dir <- args[2]
out_dir <- args[3]
plot_dir <- args[4]
k <- args[5] %>% as.integer()
z_cutoff <- args[6] %>% as.integer()
fig_width <- args[7] %>% as.integer()
fig_height <- args[8] %>% as.integer()

# chr_table_path <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/docs/chr_table.tsv"
# out_dir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/pop/lostruct/global/r"
# plot_dir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/pop/figures/lostruct/global"
# k <- 10
# z_cutoff <- 3
# fig_width <- 8
# fig_height <- 12

## Read in alternative chromosome names
chr_table <- read_tsv(chr_table_path, col_names = FALSE, col_types = cols(.default = col_character()))
if(ncol(chr_table)==1){
	chr_table <- mutate(chr_table, X2=X1)
}
chr_table <- chr_table %>%
	rename(chr=X1, lg=X2)
chrs <- chr_table$chr

##############################################################

## First, plot lostruct when all chromosomes were combined

## Read in MDS table
mds_wide <- read_tsv(str_c(in_dir, "/combined.mds.tsv"))
mds_long <- mds_wide %>%
  left_join(chr_table, by="chr") %>%
  pivot_longer(cols = str_c("mds_", 1:k), names_to = "axis", values_to = "value") %>%
  group_by(axis) %>%
  mutate(z=(value-mean(value))/sd(value)) %>%
  ungroup() %>%
  mutate(sign=case_when(z > z_cutoff ~ "plus",
                           z < -z_cutoff ~ "minus",
                           TRUE ~ "non-outlier"),
         outlier=(abs(z) > z_cutoff))
## Plot MDS
mds_plot <- mds_long %>%
	mutate(axis=fct_relevel(axis, str_c("mds_", 1:k))) %>%
	ggplot(aes(x=center/10^6, y=z, color=outlier, size=outlier)) +
	geom_point() +
  scale_color_manual(values = c("darkgrey", "magenta")) +
  scale_size_manual(values = c(0.1, 0.5)) +
	facet_grid(axis~lg, scales="free", space="free_x") +
	theme_bw() +
	xlab("position (Mbp)") +
	ylab("z-score") +
	theme(panel.spacing = unit(0.1, "lines"),
				axis.title.x=element_text(),
				legend.position="none",
				text = element_text(size=10),
				axis.text = element_text(size=6),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank()) 
ggsave(str_c(plot_dir, '/combined.mds.png'), mds_plot, width = fig_width, height = fig_height)
## Save outlier windows
mds_long %>%
  filter(outlier==TRUE) %>%
  mutate(axis_sign=str_c(axis, ".", sign)) %>%
  dplyr::select(axis_sign, prefix, chr, start, end, center) %>%
  group_by(axis_sign) %>%
  group_walk(~ write_tsv(.x, str_c(out_dir, "/combined.", .y$axis_sign, ".tsv"), col_names = FALSE))
## Save non-outlier windows
mds_long %>%
  group_by(prefix) %>%
  filter(sum(outlier)==0) %>%
  ungroup() %>%
  distinct(prefix, chr, start, end, center) %>%
  write_tsv(str_c(out_dir, "/combined.non.outlier.tsv"), col_names = FALSE)

##############################################################

## Then, plot lostruct when each chromosome was run separately

## Read in MDS table
mds_wide <- read_tsv(str_c(in_dir, "/", chrs, ".mds.tsv"))
mds_long <- mds_wide %>%
  left_join(chr_table, by="chr") %>%
  pivot_longer(cols = str_c("mds_", 1:k), names_to = "axis", values_to = "value") %>%
  group_by(axis, chr) %>%
  mutate(z=(value-mean(value))/sd(value)) %>%
  ungroup() %>%
  mutate(sign=case_when(z > z_cutoff ~ "plus",
                        z < -z_cutoff ~ "minus",
                        TRUE ~ "non-outlier"),
         outlier=(abs(z) > z_cutoff))
## Plot MDS
mds_plot <- mds_long %>%
  mutate(axis=fct_relevel(axis, str_c("mds_", 1:k))) %>%
  ggplot(aes(x=center/10^6, y=z, color=outlier, size=outlier)) +
  geom_point() +
  scale_color_manual(values = c("darkgrey", "magenta")) +
  scale_size_manual(values = c(0.1, 0.5)) +
  facet_grid(axis~lg, scales="free", space="free_x") +
  theme_bw() +
  xlab("position (Mbp)") +
  ylab("z-score") +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.title.x=element_text(),
        legend.position="none",
        text = element_text(size=10),
        axis.text = element_text(size=6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
ggsave(str_c(plot_dir, '/separated.mds.png'), mds_plot, width = fig_width, height = fig_height)
## Save outlier windows
mds_long %>%
  filter(outlier==TRUE, axis %in% str_c("mds_", 1:k)) %>%
  mutate(axis_sign=str_c(axis, ".", sign)) %>%
  dplyr::select(axis_sign, prefix, chr, start, end, center) %>%
  group_by(axis_sign, chr) %>%
  group_walk(~ write_tsv(.x, str_c(out_dir, "/", .y$chr, ".", .y$axis_sign, ".tsv"), col_names = FALSE))


# ggplot(mds_wide, aes(x=dist_3, y=dist_4, color=lg, shape=lg)) +
# 	geom_point() +
# 	scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
# 	xlab(paste0("PCo1 (", proportion_variance[1], "%)")) +
# 	ylab(paste0("PCo2 (", proportion_variance[2], "%)")) +
# 	cowplot::theme_cowplot()
# read_tsv("/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/black-and-yellow_gopher/docs/sample_table.tsv") %>%
# 	bind_cols(tibble(x=pca_summary %>% as.matrix() %>% .[2000, 4:76], y=pca_summary%>% as.matrix() %>% .[2000, 77:149])) %>%
# 	ggplot(aes(x=x, y=y, color=common_name)) +
# 	geom_point()