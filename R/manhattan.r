manhattan <- function(data, meta_data=NULL, ymin=-20, ymax=20) {

  data <- data %>% 
     rbind(meta) %>%
     rbind(meta %>% mutate(log_pvalues=-log_pvalues))
   
  don <- data %>% 
    group_by(chromosome) %>% 
    summarise(chr_len=max(position)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%  # calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(data, ., by=c("chromosome"="chromosome")) %>%  # add this to the initial dataset
    arrange(chromosome, position) %>%
    mutate(BPcum=position+tot)
 
  axisdf <- don %>% group_by(chromosome) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
  axism  <- don %>% group_by(chromosome) %>% summarise(position=max(BPcum)) 
  
  label  <- seq(ymin, ymax, length.out=5)
  limits <- c(ymin+1, ymax-1)
  breaks <- ifelse(label==0, label, ifelse(label>0, label-2, label+2))

  p <- ggplot(don, aes(x=BPcum, y=log_pvalues)) +
    geom_point(
        data =. %>% filter(type == "META"), 
        color="lightgrey", 
        size=0.8) + 
    geom_point(
        data =. %>% filter(type != "META"), 
        aes(color=as.factor(chromosome)), 
        alpha=0.8, 
        size=0.8) +
    scale_color_manual(
       values=rep(unikn::usecol(c("pal_unikn_dark"), n=4), 22)) +
    scale_x_continuous(
       name="Chromosome", 
       label=axisdf$chromosome, 
       breaks= axisdf$center, 
       minor_breaks=axism$position) +
    scale_y_continuous(
       name = "-log10(P)", 
       expand = c(0, 0), 
       limits = limits,
       label  = label, 
       breaks = breaks) +     # remove space between plot area and x axis
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  return(p)
}
