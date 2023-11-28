#!/usr/bin/env Rscript
## This script is the second dependency of /workdir/genomic-data-analysis/scripts
## It performs eigen decomposition on a covarianc matrix outputted by pcangsd

## Load required packages
suppressWarnings(suppressMessages(library(tidyverse)))

## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
prefix <- args[1]
chr <- args[2]
pc <- args[3] %>% as.integer()
cov_dir <- args[4]
beagle_dir <- args[5]
out_dir <- args[6]

# prefix <- "genome_hic_scaffold_8.w0000000000"
# chr <- "genome_hic_scaffold_8"
# pc <- 2
# cov_dir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/black-and-yellow_gopher/lostruct/global/pcangsd"
# beagle_dir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/black-and-yellow_gopher/lostruct/global/beagle"
# out_dir <- "/global/scratch/users/nicolas931010/rockfish_popgen/korean/by_group/black-and-yellow_gopher/lostruct/global/r"

## Load required packages
suppressWarnings(suppressMessages(library(tidyverse)))

## Read covariance matrix and perform eigen decomposition
suppressWarnings(suppressMessages(
  c <-read_delim(str_c(cov_dir, "/", prefix, ".cov"), delim = " ", col_names = FALSE) %>%
    as.matrix()
))
e <- eigen(c)
e_values <- e$values
e_vectors <- as_tibble(e$vectors)
## Format the output as required by lostruct 
pca_summary <- c(sum(c^2), e_values[1:pc], unlist(e_vectors[,1:pc], use.names=FALSE)) %>%
  matrix(nrow = 1) %>%
  as_tibble()
write_tsv(pca_summary, str_c(out_dir, "/", chr, ".pca_summary.tsv"), append=T, col_names=F)

## Save the start and end position of each window
suppressWarnings(suppressMessages(
  b <- read_tsv(str_c(beagle_dir, "/", prefix, ".beagle.gz"))
))
snp_position <- c(b[[1,1]], b[[dim(b)[1],1]]) %>%
  matrix(nrow = 1) %>%
  as_tibble() %>% 
  mutate(file_prefix=prefix) %>%
  relocate(file_prefix)
write_tsv(snp_position, str_c(out_dir, "/", chr, ".snp_position.tsv"), append=T, col_names=F)