rename_taxa <- function(ps2bmod, add_ASV_2tax=F){
  require(stringr)
  if (!"ASV" %in% rank_names(ps2bmod) & add_ASV_2tax==T){
    tax_table(ps2bmod) <- tax_table(cbind(ps2bmod@tax_table,rownames(ps2bmod@tax_table)))
    colnames(ps2bmod@tax_table)[length(colnames(ps2bmod@tax_table))] <- "ASV"
  }
  taxa_names(ps2bmod) <-  paste("ASV_",str_pad(as.character(1:dim(ps2bmod@otu_table)[2]),nchar(dim(ps2bmod@otu_table)[2]),pad= "0"),sep="")
  return(ps2bmod)
}