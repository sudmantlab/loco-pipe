#!/usr/bin/env Rscript
## This script plot PCA in outlier and non-outlier windows

## Load required packages
library(tidyverse)
## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
cov_dir <- args[1]
plot_dir <- args[2]
sample_table_path <- args[3]
color_by <- args[4]
chr_table_path <- args[5]
k <- args[6] %>% as.integer()

# cov_dir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/all-groups-filtered/lostruct/global/run_pcangsd_with_lostruct_outliers/"
# plot_dir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/all-groups-filtered/figures/lostruct/global"
# sample_table_path <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/all-groups-filtered/docs/sample_table.tsv"
# color_by <- "common_name"
# k <- 10
# chr_table_path <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/docs/chr_table.tsv"

## Read in alternative chromosome names
chr_table <- read_tsv(chr_table_path, col_names = FALSE, col_types = cols(.default = col_character()))
if(ncol(chr_table)==1){
  chr_table <- mutate(chr_table, X2=X1)
}
chr_table <- chr_table %>%
  rename(chr=X1, lg=X2)
chrs <- chr_table$chr

## Plot PCA with lostruct outliers when all genomes were included in a single lostruct run
sample_table <- read_tsv(sample_table_path)
cov_files <- list.files(path = cov_dir, pattern = "^combined\\..*\\.cov$") 
pca_table <- NULL
for (cov_file in cov_files){
  c <-read_delim(str_c(cov_dir, "/", cov_file), delim = " ", col_names = FALSE) %>%
    as.matrix()
  e <- eigen(c)
  e_values <- e$values
  e_vectors <- as_tibble(e$vectors) %>%
    dplyr::select(1:2) %>%
    rename(PC1=V1, PC2=V2)
  pca_table_tmp <- sample_table %>%
    bind_cols(e_vectors) %>%
    mutate(file=cov_file, 
           var_1=round(e_values[1]/sum(e_values)*100, 1),
           var_2=round(e_values[2]/sum(e_values)*100, 1))
  pca_table <- bind_rows(pca_table, pca_table_tmp)
}
pca_table_final <- pca_table %>%
  separate(col = "file", into = c("tmp_1", "axis", "sign", "tmp_2"), sep = "\\.", remove = FALSE) %>%
  dplyr::select(-tmp_1, -tmp_2) %>%
  mutate(sign_axis=str_c(axis, ".", sign, " (", var_1, "%, ", var_2, "%)"),
         sign=fct_relevel(sign, c("outlier","plus", "minus")),
         axis=fct_relevel(axis, c("non", str_c("mds_", 1:k)))) %>%
  arrange(axis, sign) %>%
  mutate(sign_axis=as_factor(sign_axis))
plot_col <- 4
plot_row <- ceiling(length(cov_files)/plot_col)
pca_plot <- pca_table_final %>%
  ggplot(aes(x=PC1, y=PC2, color=get(color_by), shape=get(color_by))) +
  geom_point() +
  scale_shape_manual(values=rep(21:25, 10)) +
  facet_wrap(~sign_axis, ncol = plot_col, scales = "free") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        panel.border = element_rect(color = "black"),
        legend.title = element_blank())
ggsave(str_c(plot_dir, '/combined.pca.png'), pca_plot, width = plot_col*2, height = plot_row*1.5)

## Plot PCA with lostruct outliers when each genome was run separately in lostruct
pca_table <- NULL
for (i in seq_len(nrow(chr_table))){
  chr <- chr_table$chr[i]
  lg <- chr_table$lg[i]
  cov_files <- list.files(path = cov_dir, pattern = str_c("^", chr, "\\..*\\.cov$"))
  for (cov_file in cov_files){
    c <-read_delim(str_c(cov_dir, "/", cov_file), delim = " ", col_names = FALSE) %>%
      as.matrix()
    e <- eigen(c)
    e_values <- e$values
    e_vectors <- as_tibble(e$vectors) %>%
      dplyr::select(1:2) %>%
      rename(PC1=V1, PC2=V2)
    pca_table_tmp <- sample_table %>%
      bind_cols(e_vectors) %>%
      mutate(chr=chr,
             lg=lg,
             file=str_remove(cov_file, str_c("^", chr, ".")), 
             var_1=round(e_values[1]/sum(e_values)*100, 1),
             var_2=round(e_values[2]/sum(e_values)*100, 1))
    pca_table <- bind_rows(pca_table, pca_table_tmp)
  }
}
pca_table_final <- pca_table %>%
  separate(col = "file", into = c("axis", "sign", "tmp"), sep = "\\.", remove = FALSE) %>%
  dplyr::select(-tmp) %>%
  mutate(sign_axis=str_c(axis, ".", sign),
         lg=as.factor(lg),
         chr=str_c(chr, " (", var_1, "%, ", var_2, "%)"),
         sign=fct_relevel(sign, c("plus", "minus")),
         axis=fct_relevel(axis, c(str_c("mds_", 1:k)))) %>%
  arrange(axis, sign, lg) %>%
  mutate(sign_axis=as_factor(sign_axis)) 
plot_col <- 4
plot_row <- ceiling(length(unique(pca_table_final$lg))/plot_col)
plot_separated_pca <- function(input){
  input %>%
    ggplot(aes(x=PC1, y=PC2, color=get(color_by), shape=get(color_by))) +
    geom_point() +
    scale_shape_manual(values=rep(21:25, 10)) +
    facet_wrap(~lg, ncol = plot_col, scales = "free", drop=FALSE) +
    theme_bw() +
    ggtitle(unique(input$sign_axis)) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top",
          panel.border = element_rect(color = "black"),
          legend.title = element_blank())
}
separated_pca <- pca_table_final %>%
  group_split(sign_axis) %>%
  lapply(plot_separated_pca)
pdf(str_c(plot_dir, '/separated.pca.pdf'), width = plot_col*2, height = plot_row*1.5+1.5)
separated_pca
dev.off()
