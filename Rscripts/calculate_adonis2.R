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
  df <- sort(pvals,decreasing = T)
  setNames(df, c("p.value","F.model","R2"))
}
calculate_adonis_jaccard2 <- function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "jaccard", binary=T)
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  R2s <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(d ~ get(colnames(ps@sam_data)[i]), df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
      Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[1]
      R2s[colnames(ps@sam_data)[i]] <- test.out$aov.tab$R2[1]
    }, error=function(e){})
  }
  df <- cbind(pvals[pvals!="NA"],Fstat[pvals!="NA"],R2s[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"])),]
  setNames(data.frame(df), c("p.value","F.model","R2"))
}
calculate_adonis_bray <- function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "bray")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  R2s <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(d ~ get(colnames(ps@sam_data)[i]), df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
      Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[1]
      R2s[colnames(ps@sam_data)[i]] <- test.out$aov.tab$R2[1]
    }, error=function(e){})
  }
  df <- cbind(pvals[pvals!="NA"],Fstat[pvals!="NA"],R2s[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"])),]
  setNames(data.frame(df), c("p.value","F.model","R2"))
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
calculate_adonis_wunifrac2 <- function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "wunifrac")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(d ~ get(colnames(ps@sam_data)[i]), df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
      Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[1]
    }, error=function(e){})
  }
  df <- cbind(pvals[pvals!="NA"],Fstat[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"])),]
  setNames(df, c("p.value","F.model","R2"))
}
calculate_adonis_unifrac <- function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "unifrac")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(d ~ get(colnames(ps@sam_data)[i]), df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
      Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[1]
    }, error=function(e){})
  }
  df <- cbind(pvals[pvals!="NA"],Fstat[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"])),]
  setNames(df, c("p.value","F.model","R2"))
}
calculate_adonis_bray3 <- function(ps, correct){
  require('vegan')
  d = phyloseq::distance(ps, "bray")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      f <- paste0("adonis(d ~ ",correct,"*",colnames(ps@sam_data)[i],", df)")
      test.out <- eval(parse(text=f))
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[2]
      Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[2]
    }, error=function(e){})
  }
  res <- cbind(pvals[pvals!="NA"],Fstat[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"])),]
  df <- res[!is.na(res[,1]),]
  setNames(df, c("p.value","F.model","R2"))
}