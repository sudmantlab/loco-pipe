## Load packages
library(tidyverse)
library(cowplot)

## Read input arguments
args <- commandArgs(trailingOnly=TRUE)
pathway <-args[1]
fig_height <- args[2] %>% as.integer()
fig_width <- args[3] %>% as.integer()

##function to normalize
norm <- function(x) x/sum(x)

sfs<- scan(paste0(outbase, '.sfs'))
p1<-barplot(sfs[-c(1,length(sfs))])
ggsave(filename = pathway, 
       plot = p1, width = fig_width, height = fig_height, units = 'in')