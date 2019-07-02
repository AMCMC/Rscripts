taxa_facet_barplot_asv <-  function(ps.temp2, Group, facet1=NULL, facet2=NULL, rank_names="ASV", lumpNA=T, taxgrp="Phylum", N=21){
  require(phyloseq)
  require(ggplot2)
  
  # Replace NA with higher level Unclassied ranks
  if (lumpNA){
    for (rank in 2:dim(ps.temp@tax_table)[2]){
      if (sum(is.na(ps.temp@tax_table[, rank])) > 0) {
        ps.temp@tax_table[is.na(ps.temp@tax_table[, rank]), rank] <-
          paste(ps.temp@tax_table[is.na(ps.temp@tax_table[, rank]), rank-1], "_Unclassified", sep = "")
      }
    }
    # Remove repetitive Unclassified.
    ps.temp@tax_table <-
      gsub("_Unclassified_.*Unclassified",
           "_Unclassified",
           ps.temp@tax_table)
  }
  
  # agglomerate on selected rank
  if (rank_names != "ASV"){ps.temp <- tax_glom(ps.temp, taxrank = rank_names, NArm = F)}


  
  # make relative & prune taxa
  ps.temp <- transform_sample_counts(ps.temp, function(x) x/sum(x))
  ps.temp <- prune_taxa(taxa_names(ps.temp) %in% names(sort(taxa_sums(ps.temp), decreasing = TRUE)[1:N]), ps.temp)
  
  # generate dataframe
  
  # relabel taxa_names
  if (rank_names=="ASV" | rank_names=="species"){
    if (c("Species") %in% rank_names(ps.temp)){
      taxa_names(ps.temp) <- make.unique(paste0(ps.temp@tax_table[,"Genus"],":",ps.temp@tax_table[,"Species"]))
    } else {taxa_names(ps.temp) <- make.unique(ps.temp@tax_table[,"Genus"])}
  } else {taxa_names(ps.temp) <- make.unique(ps.temp@tax_table[,rank_names])}
  
  # generate color vector
  bla <- as.matrix(table(ps.temp@tax_table[,taxgrp]))
  bla <- cbind(rownames(bla),bla)
  colors <- fantaxtic::gen_palette(bla)
  colvec <- unlist(colors)
  names(colvec) <- taxa_names(ps.temp)[order(paste(ps.temp@tax_table[,5], 1*taxa_sums(ps.temp)))]
  
  df <- psmelt(ps.temp)

  df$Taxa <- df$OTU
  
  # order Taxa by global abundance
  x <- aggregate(df$Abundance, by=list(df$Taxa), FUN=sum)
  #df$Taxa <- factor(df$Taxa, levels=x$Group.1[order(x$x, decreasing = F)])
  df$Taxa <- factor(df$Taxa, levels=names(colvec))
  
  # generate plot
  p <- ggplot(df, aes_string(Group, "Abundance", fill = "Taxa"))
  p <- p + geom_bar(stat = "identity")
  p <- p +  scale_fill_manual(values = colvec)
  if (!is.null(facet1) & !is.null(facet2)){
    p <- p + facet_grid(reformulate(facet1,facet2), scales="free_x", space="free")
  }
  if (!is.null(facet1) & is.null(facet2)){
    p <- p + facet_grid(reformulate(facet1), scales="free_x", space="free")
  }
  if (is.null(facet1) & !is.null(facet2)){
    p <- p + facet_grid(reformulate(".",facet2), scales="free_x", space="free")
  }

  p <- p + theme_bw() + ylab("Proportions")
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