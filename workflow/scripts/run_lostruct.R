#!/usr/bin/env Rscript
## This script read in PCA result in all windows and runs lostruct

## Load required packages
library(tidyverse)
library(lostruct)

## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
chr_table <- args[1]
pc <- args[2] %>% as.integer()
k <- args[3] %>% as.integer()
in_dir <- args[4]
out_dir <- args[5]
threads <- args[6] %>% as.integer()

# chr_list <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/docs/chr_list.txt"
# pc <- 2
# k <- 10
# out_dir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/black-and-yellow_gopher/lostruct/global/r"
# threads <- 4

## First, run lostruct with all chromosomes combined

## Read in summary tables of window-based PCA 
chrs <- read_tsv(chr_table, col_names = FALSE) %>% pull(1)
pca_summary <- read_tsv(str_c(in_dir, "/", chrs, ".pca_summary.tsv"), col_names = F) %>%
  # .[c(1:1000*10),] %>%
  as.matrix()
attr(pca_summary, 'npc') <- pc
snp_position <- read_tsv(str_c(in_dir, "/", chrs, ".snp_position.tsv"), col_names = F) %>%
  # .[c(1:1000*10),] %>%
  transmute(prefix=X1,
            chr = str_extract(X2, ".*(?=_[^_]*$)"), 
            start = as.integer(str_extract(X2, "(?<=_)[^_]*$")), 
            end = as.integer(str_extract(X3, "(?<=_)[^_]*$")), 
            center = (start+end)/2)
## Run lostruct
dist <- pc_dist(pca_summary, mc.cores = threads)
## Perform MDS with the resulting distance matrix
mds <- cmdscale(as.dist(dist), k=k, eig = T)
proportion_variance <- round(mds$eig/sum(mds$eig)*100,2)
mds_wide <- mds$points %>%
  as_tibble() %>%
  set_names(str_c("mds_", 1:k)) %>%
  bind_cols(snp_position, .)
## Write the results into text files
dist %>% 
	as_tibble() %>%
	write_tsv(str_c(out_dir, "/combined.dist.tsv"), col_names = FALSE)
write_tsv(mds_wide, str_c(out_dir, "/combined.mds.tsv"))
write_lines(proportion_variance, str_c(out_dir, "/combined.proportion_variance.txt"))

## Then, run lostruct with each chromosome separately

for (chr in chrs){
  ## Read in summary tables of window-based PCA 
  pca_summary <- read_tsv(str_c(in_dir, "/", chr, ".pca_summary.tsv"), col_names = F) %>%
    # .[c(1:1000*10),] %>%
    as.matrix()
  attr(pca_summary, 'npc') <- pc
  snp_position <- read_tsv(str_c(in_dir, "/", chr, ".snp_position.tsv"), col_names = F) %>%
    # .[c(1:1000*10),] %>%
    transmute(prefix=X1,
              chr = str_extract(X2, ".*(?=_[^_]*$)"), 
              start = as.integer(str_extract(X2, "(?<=_)[^_]*$")), 
              end = as.integer(str_extract(X3, "(?<=_)[^_]*$")), 
              center = (start+end)/2)
  ## Run lostruct
  dist <- pc_dist(pca_summary, mc.cores = threads)
  ## Perform MDS with the resulting distance matrix
  mds <- cmdscale(as.dist(dist), k=k, eig = T)
  proportion_variance <- round(mds$eig/sum(mds$eig)*100,2)
  mds_wide <- mds$points %>%
    as_tibble() %>%
    set_names(str_c("mds_", 1:k)) %>%
    bind_cols(snp_position, .)
  ## Write the results into text files
  write_tsv(mds_wide, str_c(out_dir, "/", chr, ".mds.tsv"))
  write_lines(proportion_variance, str_c(out_dir, "/", chr, ".proportion_variance.txt"))
}

