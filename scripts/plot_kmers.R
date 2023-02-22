library(ggplot2)

plot_n_kmers <- function(kmer_pheno_tbl, plot_file){


  tbl <- read.csv(kmer_pheno_tbl)


  tbl <- tbl[order(-tbl$n_kmers),]

  tbl$accession_id <- factor(tbl$accession_id, levels=tbl$accession_id)

  p <- ggplot(tbl, aes(x=accession_id, y=n_kmers, fill=as.factor(phenotype_value), label = tbl$accession_id)) + 
  
    geom_col() +
  
    annotate(geom="text",x=100, y =(max(tbl$n_kmers)/2) , label = paste("Top 3 accessions:\n1. ",tbl$accession_id[1], "\n2. ", tbl$accession_id[2], "\n3. ", tbl$accession_id[3])) +
  
    guides(fill=guide_legend(title="Phenotype value")) +
  
    theme_classic() + 
  
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
  
    scale_y_continuous(expand = c(0, 0)) +
    
    ylab("number of kmers")

  
  ggsave(plot_file, p)
  
}

plot_n_kmers(snakemake@input[[1]], snakemake@output[[1]])
       
