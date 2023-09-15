plot_fst_per_snp <- function(x){
  ggplot(x, aes(x=pos/10^6, y=fst)) +
    geom_point(aes(color = lg), size=0.1, alpha=0.5) +
    geom_smooth(color='blue', se=F) +
    scale_color_manual(values = rep(c('black', 'grey50'), 50)) +
    scale_x_continuous(breaks=seq(0, 100, 10)) +
    coord_cartesian(ylim=c(0, 1), expand = F) +
    xlab('position (Mbp)') +
    ylab('Fst') +
    facet_grid(.~lg, scales='free_x', space='free_x') +
    theme_cowplot() +
    theme(panel.spacing = unit(0.1, 'lines'),
          axis.title.x=element_text(),
          legend.position='none',
          text = element_text(size=10),
          axis.text = element_text(size=6))
}
fixed_snp_window_fst <- function(x, window_length){
  x %>%
    group_by(lg) %>%
    mutate(index = row_number(), 
           index = ceiling(index / window_length)) %>%
    group_by(lg, index) %>%
    filter(n() == window_length) %>%
    summarise(pos = (min(pos)+max(pos))/2, fst=sum(alpha)/sum(beta)) %>%
    ungroup()
}
fixed_length_window_fst <- function(x, window_length){
  x %>%
    mutate(pos=cut(pos, breaks=seq(0,100*10^6,window_length), labels=seq(window_length/2,100*10^6-window_length/2,window_length))) %>%
    group_by(lg, pos) %>%
    summarise(fst=sum(alpha)/sum(beta), n_snps=n()) %>%
    mutate(pos=as.numeric(as.character(pos)))
}