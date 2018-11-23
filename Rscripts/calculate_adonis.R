calculate_adonis_jaccard <- function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "jaccard", binary=T)
  df = as(phyloseq::sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(d ~ get(colnames(ps@sam_data)[i]), df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
    }, error=function(e){})
  }
  names(pvals) <- colnames(ps@sam_data)
  sort(pvals,decreasing = T)
}
calculate_adonis_bray <- function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "bray")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(d ~ get(colnames(ps@sam_data)[i]), df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
    }, error=function(e){})
  }
  names(pvals) <- colnames(ps@sam_data)
  sort(pvals,decreasing = T)
}
calculate_adonis_wunifrac <- function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "wunifrac")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(d ~ get(colnames(ps@sam_data)[i]), df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
    }, error=function(e){})
  }
  names(pvals) <- colnames(ps@sam_data)
  sort(pvals,decreasing = T)
}