#### get all node descendant sums and corresponding consensus taxonomy given a phyloseq object with a tree ####
# preferably rarefy ps object before

ps.temp <- ps

#if the tree is not rooted mid point root the tree
if (!ape::is.rooted(ps.temp@phy_tree)){
  phy_tree(ps.temp) <- phytools::midpoint.root(ps.temp@phy_tree)}

#for each node get all tips
node.desc <- phangorn::Descendants(ps.temp@phy_tree, node = (1:ps.temp@phy_tree$Nnode)+ntaxa(ps.temp), "tips")
# give an arbitrary label
names(node.desc) <- paste0("Node:", 1:ps.temp@phy_tree$Nnode)
# get number of descendent tips
node.desc.num_tips <- lapply(node.desc, function(x) length(x))
# for each node get the sum off all descendants for each sample

# get the associated tip labels
node.desc.labels <- lapply(node.desc, function(x) ps.temp@phy_tree$tip.label[x])

# get the node sum matrix
node.desc.sums <- lapply(node.desc.labels, function(x) rowSums(ps.temp@otu_table[,unlist(x)]))
node.desc.sums <- plyr::ldply(node.desc.sums)
rownames(node.desc.sums) <- node.desc.sums[,1]
node.desc.sums <- node.desc.sums[,-1]

#get node consensus taxonomy (doesnt handle complete non uniques; ei root)
consensus_tax <- function(lot){
  tryCatch( paste(unique(lot[,apply(lot, 2, function(x) length(unique(x)))==1]), collapse = ";"),
            error=function(e) NA)
}

node.desc.tax <- unlist(lapply(node.desc.labels, function(x) consensus_tax(ps@tax_table[x,])))