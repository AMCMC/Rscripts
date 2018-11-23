buildps <- function(project, includemeta="", buildtree=F){
  
#  project <- "2018_08_Preterm_Sepsis"
#  additional_meta <- ""
#  buildtree=F

  meta_data <- read.csv("C:/Users/mdavi/Documents/AMC/MicrobiotaCentre/Libraries/Overview metadata seq.txt",sep="\t", header=T, skip=4)
  meta_data <- meta_data[!meta_data$Seq_ID %in% c("","N/A"),]
  table(droplevels(meta_data[meta_data$Project_name==project,]$Seq_ID))
  
  ps <- list()
  ps.p <- list()
  for (i in names(table(droplevels(meta_data[meta_data$Project_name==project,]$Seq_ID)))){
    ps[[i]] <- paste0("../../Libraries/ps_objects/",i,"clean.RDS")
    ps.p[[i]] <- readRDS(paste0("../../Libraries/ps_objects/",i,"clean.RDS"))
    ps.p[[i]] <- prune_samples(ps.p[[i]]@sam_data$Project_name==project,ps.p[[i]])
    ps.p[[i]] <- prune_taxa(names(which(taxa_sums(ps.p[[i]])!=0)),ps.p[[i]])
    sample_names(ps.p[[i]]) <- ps.p[[i]]@sam_data$Int_ID
    sample_names(ps.p[[i]]) <- ps.p[[i]]@sam_data$Ext_ID
  }
  
  
  for (p in names(ps.p)){
    ps.p[[p]]@phy_tree <- NULL
    taxa_names(ps.p[[p]]) <- ps.p[[p]]@tax_table[,"ASV"]
  }
  
  ps.master <- list()
  if (length(ps.p)==1){
    ps.master[["ps.sub"]] <- ps.p[[1]]
  } else{
    for (i in 1:length(ps.p)){
      ps.p[[i]]@phy_tree <- NULL
      taxa_names(ps.p[[i]]) <- ps.p[[i]]@tax_table[,"ASV"]
    }
    ps.master[["ps.sub"]] <- merge_phyloseq(ps.p[[1]], ps.p[[2]])
    if (length(ps.p)>2){
      for (i in 3:length(ps.p)){
        ps.master[["ps.sub"]] <- merge_phyloseq(ps.master[["ps.sub"]], ps.p[[i]])
      }
    }
  }
  
  taxa_names(ps.master[["ps.sub"]]) <-  paste("ASV_",str_pad(as.character(1:dim(ps.master[["ps.sub"]]@otu_table)[2]),nchar(dim(ps.master[["ps.sub"]]@otu_table)[2]),pad= "0"),sep="")
  
  if (additional_meta!=""){
    meta <- read.csv(additional_meta,sep="\t", row.names = 1)
    sample_data(ps.master[["ps.sub"]]) <- sample_data(cbind(sample_data(ps.master[["ps.sub"]]),meta[ps.master[["ps.sub"]]@sam_data$Ext_ID,]))
  }
  
  
  #ps@phy_tree <- phy_tree(build_tree_from_asv(ps))
  #fitGTR <- readRDS("dada2.fitGTR.RDS")
  if (buildtree==T){
    phy_tree(ps.master[["ps.sub"]]) <- phy_tree(build_tree_from_asv(ps.master[["ps.sub"]], optim = F))
  }
  
  return(ps)
}
