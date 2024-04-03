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
## getting the chromosomes
chr<-read_lines(chr_list,skip_empty_rows = TRUE)


## testing
#indir <-"/global/scratch/users/zzhou32/softwares/loco-pipe/toyfish/angsd/get_theta"
#outdir <- "/global/scratch/users/zzhou32/softwares/loco-pipe/toyfish/figures/theta/sfs_distribution"
#population<- "sunset"
#chr_list<- "/global/scratch/users/zzhou32/softwares/loco-pipe/toyfish/docs/chr_list.txt"
#ref_type<-1
#fig_height <- args[5] %>% as.integer()
#fig_width <- args[6] %>% as.integer()

## getting the dataframe, and slicing in half if it is folded
sfs_list<-read_delim(str_c(indir,"/",population,".",chr,".sfs"), col_names=FALSE, delim=" ")
n<- (ncol(sfs_list)-2)/2
## calculating the expected sfs in a neutrally evolving population
proportion_vector<-1/sum(1/(1:(2*n-1))) / (1:(2*n-1))
summed<- proportion_vector+rev(proportion_vector)

## transforming the sfs_list into a dataframe
df <-sfs_list %>% 
  select(if (ref_type == 1) 1:(floor(ncol(.)+1)/ 2) else 1:(floor(ncol(.))-1)) %>% 
  summarise(across(everything(), ~sum(., na.rm = TRUE))) %>% 
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(sf = str_extract(Variable, "\\d+") %>% as.integer()) %>%
  mutate(sf=sf-1) 

if (ref_type == 1){
  df_trim <- df %>% slice(-1)
  foldness <- "Number of Minor Alleles"
  hypothetical_distribution<-summed[1:n]
  hypothetical_distribution[[n]]=hypothetical_distribution[[n]]/2
} else{
  df_trim <-df %>% slice(-c(1, nrow(.)))
  foldness <- "Number of Derived Alleles"
  hypothetical_distribution<-proportion_vector
}

## for variable_sites_only, we calculate the proportion
df_trim <- df_trim %>% mutate(estimated_proportion= hypothetical_distribution) %>%
  mutate(real_proportion=Value/sum(Value))

## managing the position of the annotations
max_y<-max(max(df_trim$estimated_proportion),max(df_trim$real_proportion))
max_x<-max(df_trim$sf)

with_invariable<-ggplot(df, aes(x=sf, y=Value)) + geom_bar(stat="identity", fill="darkblue")+
  ggtitle("including invariable sites")+xlab(foldness)+ 
  ylab("Counts (in log scale)")+scale_y_log10()+theme_bw()

variable_only<-ggplot(df_trim, aes(x=sf))+geom_bar(aes(y=real_proportion),stat="identity",fill="darkred")+
  geom_line(aes(y=estimated_proportion,group=1),color="blue")+
  geom_point(aes(y=estimated_proportion),color="blue")+
  ggtitle("variable sites only")+xlab(foldness)+ylab("Proportion")+
  annotate("text",label="empirical SFS estimated from ANGSD",color="darkred",x=max_x*5/9,y=max_y*7/9)+
  annotate("text",label="expected SFS in a neutrally evolving population",color="blue",x=max_x*5/9,y=max_y*6/9)+
  theme_bw()
variable_only


p1<-plot_grid(with_invariable,variable_only)
 
p1

ggsave(filename = paste0(outdir,'/', population, '.','sfs.png'), 
       plot = p1, width = 10, height = 5, units = 'in')
