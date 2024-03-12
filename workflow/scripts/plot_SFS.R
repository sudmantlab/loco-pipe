## Load packages
library(tidyverse)
library(cowplot)
library(ggplot2)

## Read input arguments
args <- commandArgs(trailingOnly=TRUE)
indir <-args[1]
outdir <-args[2]
population<-args[3]
chr_list<-args[4]
ref_type<-args[5] %>% as.integer()

## testing
#indir <-"/global/scratch/users/zzhou32/softwares/loco-pipe/toyfish/angsd/get_theta"
#outdir <- "/global/scratch/users/zzhou32/softwares/loco-pipe/toyfish/figures/theta/sfs_distribution"
#population<- "sunset"
#chr_list<- "/global/scratch/users/zzhou32/softwares/loco-pipe/toyfish/docs/chr_list.txt"
ref_type<-1
#fig_height <- args[5] %>% as.integer()
#fig_width <- args[6] %>% as.integer()

## getting the chromosomes
chr<-read_lines(chr_list,skip_empty_rows = TRUE)

## getting the dataframe, and slicing in half if it is folded
df <-read_delim(str_c(indir,"/",population,".",chr,".sfs"), col_names=FALSE, delim=" ") %>% 
  select(if (ref_type == 1) 1:(floor(ncol(.)+1)/ 2) else everything()) %>% 
  summarise(across(everything(), ~sum(., na.rm = TRUE))) %>% 
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(sf = str_extract(Variable, "\\d+") %>% as.integer()) %>%
  mutate(sf=sf-1)

if (ref_type == 1){
  df_trim <- df %>% slice(-1)
  foldness <- "Number of Minor Alleles"
} else{
  df_trim <-df %>% slice(-c(1, nrow(.)))
  foldness <- "Number of Derived Alleles"
}


with_invariable<-ggplot(df, aes(x=sf, y=Value)) + geom_bar(stat="identity", fill="darkblue")+
  ggtitle("including invariable sites")+xlab(foldness)+ 
  ylab("Counts (in log scale)")+scale_y_log10()

variable_only<-ggplot(df_trim, aes(x=sf, y=Value))+geom_bar(stat="identity",fill="darkred")+
  ggtitle("variable sites only")+xlab(foldness)+ylab("Counts")

p1<-plot_grid(with_invariable,variable_only)
p1
ggsave(filename = paste0(outdir,'/', population, '.','sfs.png'), 
       plot = p1, width = 10, height = 5, units = 'in')