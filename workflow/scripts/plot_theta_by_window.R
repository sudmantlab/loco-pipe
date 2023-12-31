## Load packages
library(tidyverse)
library(cowplot)

## Read input arguments
args <- commandArgs(trailingOnly=TRUE)
indir <- args[1]
plot <-  args[2]
chr_table <- args[3]
window_size <- args[4] %>% as.integer()
step_size <- args[5] %>% as.integer()
population <- args[6]
fig_height <- args[7] %>% as.integer()
fig_width <- args[8] %>% as.integer()

#indir <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/angsd/get_theta"
#plot <-  "/global/scratch/users/nicolas931010/loco-pipe/toyfish/figures/get_theta/vermilion.theta_by_window.png"
#chr_table <- "/global/scratch/users/nicolas931010/loco-pipe/toyfish/docs/chr_table.tsv"
#window_size <- 10000
#step_size <- 10000
#population <- "vermilion"

## List of chromosomes
chr_table <- read_tsv(chr_table, col_names = FALSE, col_types = cols(.default = col_character()))
if (ncol(chr_table) == 1){
  chr_table <- mutate(chr_table, X2 = X1)
}
chr_table <- chr_table %>%
  rename(chr = X1, lg = X2) #rename the columns for the plot
chrs <- chr_table$chr
## read in the theta result
theta_table <- read_tsv(paste0(indir, '/', population, '.', chrs, '.', window_size, 'window_', step_size, 'step.thetas.pestPG')) 
theta_table_long <- theta_table %>%
  left_join(chr_table, by=c("Chr"="chr")) %>%
  transmute(lg=lg, pos=WinCenter, pi=tP/nSites, `Watterson's theta`=tW/nSites, `Tajima's D`=Tajima) %>%
  pivot_longer(cols = 3:5, names_to = "name", values_to = "value") %>%
  mutate(name=fct_relevel(name, c("pi", "Watterson's theta" , "Tajima's D")))
## plot theta in sliding windows
theta_plot <- theta_table_long %>%
  ggplot(aes(x=pos/10^6, y=value)) +
  geom_point(size = 0.1) +
  facet_grid(name~lg, space = "free_x", scale="free") +
  xlab('position (Mbp)') +
  theme_cowplot() +
  theme(panel.spacing = unit(0.1, 'lines'),
        axis.title.x = element_text(),
        legend.position = 'none',
        text = element_text(size = 10),
        axis.text = element_text(size = 6), 
        panel.border = element_rect(color="black"),
        axis.line = element_blank())
ggsave(filename = plot, 
       plot = theta_plot, width = fig_width, height = fig_height, units = 'in', pointsize = 20)
