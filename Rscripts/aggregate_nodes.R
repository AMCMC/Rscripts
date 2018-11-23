ps_node_descendants <- function(ps){
  
  for (node in (ps.temp@phy_tree$Nnode+2):max(ps.temp@phy_tree$edge)){
    node.desc[[paste0("X",node)]] <- unlist(phangorn::Descendants(ps.temp@phy_tree, node = node, "tips"))
    #node.anc[[node]] <- phangorn::Ancestors(ps.temp@phy_tree, node = node, "all")
  }
}



names(which(table(ps@phy_tree$edge[,1])!=1))


ps_node_count_aggregation <- function(ps){
  require(ape)
  require(phytools)
  
  node.desc <- list()
  node.anc <- list()
  
  ps.temp <- ps
  if (!ape::is.rooted(ps.temp@phy_tree)){
  phy_tree(ps.temp) <- phytools::midpoint.root(ps.temp@phy_tree)
  } else {
    
  }
  
  for (node in (ps.temp@phy_tree$Nnode+2):max(ps.temp@phy_tree$edge)){
    node.desc[[paste0("X",node)]] <- unlist(phangorn::Descendants(ps.temp@phy_tree, node = node, "tips"))
    #node.anc[[node]] <- phangorn::Ancestors(ps.temp@phy_tree, node = node, "all")
  }
  
  node_sums <- lapply(node.desc, function(x) rowSums(ps.temp@otu_table[,unlist(x)]))
  node_sums <- plyr::ldply(node_sums)
  rownames(node_sums) <- node_sums[,1]
  node_sums <- node_sums[,-1]
  
  return(node_sums)
}

ps_nodeconsensus_taxonomy <- function(ps){
  require(ape)
  require(phytools)
  
  node.desc <- list()
  node.anc <- list()
  
  ps.temp <- ps
  if (!ape::is.rooted(ps.temp@phy_tree)){
    phy_tree(ps.temp) <- phytools::midpoint.root(ps.temp@phy_tree)
  } else {
    
  }
  
  for (node in (ps.temp@phy_tree$Nnode+2):max(ps.temp@phy_tree$edge)){
    node.desc[[paste0("X",node)]] <- unlist(phangorn::Descendants(ps.temp@phy_tree, node = node, "tips"))
    #node.anc[[node]] <- phangorn::Ancestors(ps.temp@phy_tree, node = node, "all")
  }
  
  node_sums <- lapply(node.desc, function(x) rowSums(ps.temp@otu_table[,unlist(x)]))
  node_sums <- plyr::ldply(node_sums)
  rownames(node_sums) <- node_sums[,1]
  node_sums <- node_sums[,-1]
  
  return(node_sums)
}