taxa_facet_barplot_asv <-  function(ps, Group="Subject_ID", facet1=NULL, facet2=NULL, taxrank="ASV", lumpNA=T, topglobal=T, relative=T, collapse_species=10){
  #ps; Group="Subject_ID"; facet1=NULL; facet2=NULL; taxrank="ASV"; lumpNA=T; topglobal=T; relative=T;collapse_species=2
  require(phyloseq)
  require(ggplot2)
  if (collapse_species<999){ps@tax_table[,"Species"][lapply(strsplit(ps@tax_table[,7],"/"), length)>collapse_species] <- "Ambigous"}
  
  # Replace NA with higher level Unclassied ranks
  if (lumpNA){
    for (rank in 2:dim(ps@tax_table)[2]){
      if (sum(is.na(ps@tax_table[, rank])) > 0 & rank_names(ps)[rank]!="Species") {
        ps@tax_table[is.na(ps@tax_table[, rank]), rank] <-
          paste(ps@tax_table[is.na(ps@tax_table[, rank]), rank-1], "_Unclassified", sep = "")
      }
    }
    # Remove repetitive Unclassified.
    ps@tax_table <-
      gsub("_Unclassified_.*Unclassified",
           "_Unclassified",
           ps@tax_table)
  }
  
  # color vector to use
  colours <-
    c(
      "#F0A3FF","#0075DC","#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5",
      "#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010",
      "#5EF1F2","#00998F","#740AFF","#990000","#FFFF00"
    )
  
  # agglomerate on selected rank
  if (taxrank != "ASV"){ps <- tax_glom(ps, taxrank = taxrank, NArm = F)}
  
  # How many clades to plot
  N <- min(length(taxa_names(ps)), length(colours))
  
  # make relative & prune taxa
  if (relative){ps <- transform_sample_counts(ps, function(x) x/sum(x))}
  if (topglobal){
    selecttaxa <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:N])
    ps <- prune_taxa(taxa_names(ps) %in% selecttaxa, ps)
  } else {
    selecttaxa <- taxa_names(ps)[order(colMaxs(ps@otu_table), decreasing = T)[1:N]]
    ps <- prune_taxa(taxa_names(ps) %in% selecttaxa, ps)
  }
  
  # generate dataframe
  df <- psmelt(ps)
  df$Taxa <- df[,taxrank]
  if (taxrank=="ASV" | taxrank=="Species"){if (c("Species") %in% rank_names(ps)){
    df$Taxa <- paste0(df$OTU,":",df$Genus,":",df$Species)} else {
    df$Taxa <- paste0(df$OTU,":",df$Genus)}}
  
  # order Taxa by global abundance
  x <- aggregate(df$Abundance, by=list(df$Taxa), FUN=sum)
  df$Taxa <- factor(df$Taxa, levels=x$Group.1[order(x$x, decreasing = F)])
  
  # generate plot
  p <- ggplot(df, aes_string(Group, "Abundance", fill = "Taxa"))
  p <- p + geom_bar(stat = "identity")
  p <- p +  scale_fill_manual(values = c(rev(colours[0:(N)])))
  
  if (!is.null(facet1) & !is.null(facet2)){
    p <- p + facet_grid(reformulate(facet1,facet2), scales="free_x", space="free")
  }
  if (!is.null(facet1) & is.null(facet2)){
    p <- p + facet_grid(reformulate(facet1), scales="free_x", space="free")
  }
  if (is.null(facet1) & !is.null(facet2)){
    p <- p + facet_grid(reformulate(".",facet2), scales="free_x", space="free")
  }
  
  p <- p + theme_bw()
  if (relative){p <- p + ylab("Proportions")}else{p <- p + ylab("Counts")}
  
  p <- p + scale_y_continuous(expand = c(0, 0)) +
    theme(strip.background = element_rect(fill = "gray85")) +
    theme(panel.spacing = unit(0.3, "lines"))
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p <- p + guides(fill = guide_legend(ncol = 1, reverse = F))
  p <- p + theme(panel.spacing = unit(1, "lines"))
  p <- p + theme(text = element_text(size = 20),
                 axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}
