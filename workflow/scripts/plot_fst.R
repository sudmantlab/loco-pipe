## Load packages
library(tidyverse)
library(cowplot)

## Read input arguments
args <- commandArgs(trailingOnly=TRUE)
chr_table <- args[1]
basedir <- args[2]
population1 <- args[3]
population2 <- args[4]
snp_window_size <- args[5] %>% as.integer()
bp_window_size <- args[6] %>% as.integer()
fig_height <- args[7] %>% as.integer()
fig_width <- args[8] %>% as.integer()

## List of chromosomes
chr_table <- read_tsv(chr_table, col_names = FALSE, col_types = cols(.default = col_character()))
if (ncol(chr_table) == 1){
  chr_table <- mutate(chr_table, X2 = X1)
}
chr_table <- chr_table %>%
  rename(chr = X1, lg = X2) #rename the columns for the plot
chrs <- chr_table$chr

## Define some useful functions

# This functions makes the Fst plot with SNP-windows.
plot_fst <- function(x){
  ggplot(x, aes(x = pos/10^6, y = fst)) +
    geom_point(size = 0.1, alpha = 0.5) +
    ## uncomment the following to alternate colors across chromosomes
    #geom_point(aes(color = lg), size = 0.1, alpha = 0.5) +
    #scale_color_manual(values = rep(c('black', 'grey50'), 50)) +
    geom_smooth(color = 'blue', se = F) +
    ## uncomment the following to manually set x-axis breaks
    #scale_x_continuous(breaks = seq(0, 100, 10)) +
    coord_cartesian(ylim = c(0, 1), expand = F) +
    xlab('position (Mbp)') +
    ylab('Fst') +
    facet_grid(.~lg, scales = 'free_x', space = 'free_x') +
    theme_cowplot() +
    theme(panel.spacing = unit(0.1, 'lines'),
          axis.title.x = element_text(),
          legend.position = 'none',
          text = element_text(size = 10),
          axis.text = element_text(size = 6), 
          panel.border = element_rect(color="black"),
          axis.line = element_blank())
}

# This function subsets each linkage group via SNP window. (The window_length is the number of 
# SNPs contained in one window.) 
fixed_snp_window_fst <- function(x, window_length){
  x %>%
    group_by(lg) %>%
    mutate(index = row_number(), 
           index = ceiling(index / window_length)) %>%
    group_by(lg, index) %>%
    filter(n() == window_length) %>%
    summarise(start_pos = min(pos), end_pos = max(pos), alpha = sum(alpha), beta = sum(beta)) %>%
    ungroup() %>%
    mutate(pos = (start_pos+end_pos)/2, fst = alpha/beta)
}

# This function subsets each linkage group via base pair (BP) window. (The window_length is the number of 
# base pairs contained in one window.) 
fixed_bp_window_fst <- function(x, window_length){
  x %>%
    mutate(pos = cut(pos, breaks = seq(0, 100*10^6, window_length), 
                     labels = seq( window_length/2, 100*10^6-window_length/2, window_length))) %>%
    group_by(lg, pos) %>%
    summarise(alpha = sum(alpha), beta = sum(beta), n_snps = n()) %>%
    ungroup() %>%
    mutate(fst = alpha/beta,
           pos = as.numeric(as.character(pos)),
           start_pos = pos-window_length/2,
           end_pos = pos+window_length/2)
}
## read in Fst
fst <- read_tsv(paste0(basedir, '/angsd/get_fst/', population1, '.', population2, '.', chrs, '.alpha_beta.txt'), 
                col_names = c('chr', 'pos', 'alpha', 'beta')) %>%
  mutate(fst = alpha/beta) %>%
  left_join(chr_table, by = "chr")

## Plot Fst in fixed SNP windows
fst_snp_window <- fixed_snp_window_fst(fst, snp_window_size)
write_tsv(fst_snp_window, paste0(basedir, '/angsd/get_fst/', population1, '.', population2, '.', snp_window_size, 'snp_window.tsv'))
fst_snp_window_plot <- plot_fst(fst_snp_window)
ggsave(filename  = paste0(basedir, '/figures/fst/', population1, '.', population2, '.', snp_window_size, 'snp_window.png'), 
       plot = fst_snp_window_plot, width = fig_width, height = fig_height, units = 'in', pointsize = 20)

## Plot Fst in fixed bp windows
fst_bp_window <- fixed_bp_window_fst(fst, bp_window_size)
write_tsv(fst_bp_window, paste0(basedir, '/angsd/get_fst/', population1, '.', population2, '.', bp_window_size, 'bp_window.tsv'))
fst_bp_window_plot <- plot_fst(fst_bp_window)
ggsave(filename  = paste0(basedir, '/figures/fst/', population1, '.', population2, '.', bp_window_size, 'bp_window.png'), 
       plot= fst_bp_window_plot, width = fig_width, height = fig_height, units = 'in', pointsize = 20)

## Plot original Fst with no windows
fst_plot <- plot_fst(fst)
ggsave(filename  = paste0(basedir, '/figures/fst/', population1, '.', population2, '.png'), 
       plot = fst_plot, width = fig_width, height = fig_height, units = 'in', pointsize = 20)
