boxplot_ASV_list <- function(asvlist, ps, x) {
  pst <- psmelt(prune_taxa(taxa_names(ps) %in% asvlist, ps))
  p <- ggplot(pst, aes_string(x=x, y="Abundance +1")) + 
    geom_boxplot() + 
    ggbeeswarm::geom_beeswarm() + 
    facet_wrap(~OTU) +
    scale_y_log10() + 
    NULL
  return(p)
}