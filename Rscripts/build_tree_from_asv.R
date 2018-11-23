build_tree_from_asv <- function(ps, optim=T){
  require("phangorn")
  require("DECIPHER")
  seqs <- as.vector(ps@tax_table[,8])
  names(seqs) <- taxa_names(ps)
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm)
  
  fit = pml(treeNJ, data=phang.align)
  fitGTR <- update(fit, k=4, inv=0.2)
  if (optim==T){
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  }
  fitGTR$tree
}