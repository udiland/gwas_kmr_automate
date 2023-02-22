
library(dplyr)
library(ggplot2)
library(htmlwidgets)
library(plotly)

make_manhattan <- function(in_pass_5, plot_all_pass5, plot_inter_pass5){

  # read results
  kmer_tbl <- read.csv(in_pass_5)
  
  # take only kmers tha were aligned
  kmers <- kmer_tbl[kmer_tbl$state != "not_aligned" ,]


  # mark partial math of alignment
  kmers$match <- ifelse(kmers$match == "full", kmers$match, "partial")

  # mark if the kmer is match uniqly
  kmers$state <- ifelse(kmers$state == "duplicate", kmers$state, "uniqe")


  # create table for plotting all kmers
  don <- kmers %>% group_by(chr_kmer) %>% summarise(chr_len=max(end_kmer)) %>%
  
    #calculate comulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len) %>%
  
    # Add this info to the initial dataset
    left_join(kmers, ., by=c("chr_kmer"="chr_kmer")) %>%
  
    # Add a cumulative position of each SNP
    arrange(chr_kmer, end_kmer) %>%
    mutate( BPcum=end_kmer+tot)

  
  # create a table with midlle position for each chromosome (for x axis plot) 
  axisdf = don %>% group_by(chr_kmer) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


  don$text <- paste("kmer: ", don$name_kmer, "\nStart: ", don$start_kmer, 
                  "\nEnd: ", don$end_kmer, "\nChromosome: ", don$chr_kmer, 
                  "\nlogPV:",don$logPV %>% round(6),
                  "\nScaffold: ", don$name_contig,
                  "\nStart: ", don$start_contig,
                  "\nEnd: ", don$end_contig,
                  sep = "")

  don$link <- paste("http://icci1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                  don$chr_kmer, "%3A" ,don$start_kmer, "%2D", don$end_kmer,
                  "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")


  p <- ggplot(don, aes(x=BPcum, y=logPV, text=text, customdata=link, shape=state)) +
  
    # Show all points  (duplicate kmers are circles and uniq are triengels)
    geom_jitter(aes(color=as.factor(match)), size=2) +
  
    # set the shapes
    scale_shape_manual(values = c(16,8)) +
  
    # custom X axis:
    scale_x_continuous( label = axisdf$chr_kmer, breaks= axisdf$center ) +
  
    # remove space between plot area and x axis and set the y axis limits
    scale_y_continuous(expand = c(0, 0), limits = c(min(don$logPV - 0.03), max(don$logPV) + 0.05)) +
  
    # Custom the theme:
    theme_bw() +
    theme(
      #legend.title=element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x=element_blank()
    ) +
  
    # add title from the user
    ggtitle("kmers Manhattan") +
    theme(plot.title = element_text(hjust = 0.5, size = 22))

  p <- p + theme(legend.title = element_blank()) 

  p2 <- ggplotly(p, tooltip = c('text'))

  p2$x$data[[1]]$customdata <- don$link

  js <- "
  function(el) {
    el.on('plotly_click', function(d) {
      var url = d.points[0].customdata;
      //url
      window.open(url);
    });
  }"


  p3 <- onRender(p2, js)

  l <- plotly::ggplotly(p3)
  
  htmlwidgets::saveWidget(l, plot_all_pass5)

  ######################################################################################

  # take only kmers that intersect with scaffold
  kmers <- kmers %>% group_by(name_kmer) %>% filter(startsWith(as.character(name_contig), "NODE"))

  # create table for plotting
  don <- kmers %>% group_by(chr_kmer) %>% summarise(chr_len=max(end_kmer)) %>%
  
    #calculate comulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len) %>%
  
    # Add this info to the initial dataset
    left_join(kmers, ., by=c("chr_kmer"="chr_kmer")) %>%
  
    # Add a cumulative position of each SNP
    arrange(chr_kmer, end_kmer) %>%
    mutate( BPcum=end_kmer+tot)

  # create a table with midlle position for each chromosome (for x axis plot)
  axisdf = don %>% group_by(chr_kmer) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


  don$text <- paste("kmer: ", don$name_kmer, "\nStart: ", don$start_kmer, 
                  "\nEnd: ", don$end_kmer, "\nChromosome: ", don$chr_kmer, 
                  "\nlogPV:",don$logPV %>% round(6),
                  "\nScaffold: ", don$name_contig,
                  "\nStart: ", don$start_contig,
                  "\nEnd: ", don$end_contig,
                  sep = "")

  don$link <- paste("http://icci1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                  don$chr_kmer, "%3A" ,don$start_kmer, "%2D", don$end_kmer,
                  "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")


  p <- ggplot(don, aes(x=BPcum, y=logPV, text=text, customdata=link, shape=state)) +
  
    # Show all points  (duplicate kmers are circles and uniq are triengels)
    geom_jitter(aes(color=as.factor(match)), size=2) +
  
    #scale_color_manual(values = rep(c("blue", "orange"), 4 )) +
  
    # set the shapes
    scale_shape_manual(values = c(16,8)) +
  
    # custom X axis:
    scale_x_continuous( label = axisdf$chr_kmer, breaks= axisdf$center ) +
  
    # remove space between plot area and x axis and set the y axis limits
    scale_y_continuous(expand = c(0, 0), limits = c(min(don$logPV - 0.03), max(don$logPV) + 0.05)) +
  
  
    # Custom the theme:
    theme_bw() +
    theme(
      #legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x=element_blank()
    ) +
  
    # add title from the user
    ggtitle("kmers Manhattan interscet only") +
    theme(plot.title = element_text(hjust = 0.5, size = 22))

  p <- p + labs(color="") 

  p2 <- ggplotly(p, tooltip = c('text'))

  p2$x$data[[1]]$customdata <- don$link

  js <- "
  function(el) {
    el.on('plotly_click', function(d) {
      var url = d.points[0].customdata;
      //url
      window.open(url);
    });
  }"


  p3 <- onRender(p2, js)

  l <- plotly::ggplotly(p3)
  
  htmlwidgets::saveWidget(l, plot_inter_pass5)
}

make_manhattan(snakemake@input[[1]], snakemake@output[[1]], snakemake@output[[2]])
