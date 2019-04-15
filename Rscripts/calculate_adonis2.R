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
calculate_adonis_bray_md <- function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "bray")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  R2s <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      var <- colnames(df)[i]
      d2 <- as.matrix(d)[!(is.na(df[,i])),!(is.na(df[,i]))]
      df2 <- df[!(is.na(df[,i])),]
      test.out <- adonis(d2 ~ get(var), df2)
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
calculate_adonis_philr <- function(ps){
  require('philr')
  require('vegan')
  # warning rooted?
  # warning branhing?
  suppressWarnings(
  ps.sub.philr <- philr(otu_table(ps.sub) + 1, phy_tree(ps.sub), 
                        part.weights='enorm.x.gm.counts', 
                        ilr.weights='blw.sqrt'))
  
  d = dist(ps.sub.philr, method="euclidean")
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
calculate_adonis_philr_drop_NA <- function(ps){
  require('philr')
  require('vegan')
  # warning rooted?
  # warning branhing?
  suppressWarnings(
    ps.sub.philr <- philr(otu_table(ps.sub) + 1, phy_tree(ps.sub), 
                          part.weights='enorm.x.gm.counts', 
                          ilr.weights='blw.sqrt'))
  
  d = dist(ps.sub.philr, method="euclidean")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  R2s <- pvals
  samples_used <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      samples_used[colnames(ps@sam_data)[i]] <- sum(!is.na(df[,i]))
      test.out <- adonis(as.matrix(d)[!is.na(df[,"Negatief"]),!is.na(df[,"Negatief"])] ~ get(colnames(ps@sam_data)[i])[!is.na(df[,"Negatief"])], df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
      Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[1]
      R2s[colnames(ps@sam_data)[i]] <- test.out$aov.tab$R2[1]
    }, error=function(e){})
  }
  df <- cbind(samples_used[pvals!="NA"],pvals[pvals!="NA"],Fstat[pvals!="NA"],R2s[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"])),]
  setNames(data.frame(df), c("samples_used","p.value","F.model","R2"))
}
calculate_adonis <- function(ps, distance_m){
  options(error=NULL)
  require('philr')
  require('vegan')
  # warning rooted?
  # warning branhing?
  # warning normalize
  supported_methods <- c(c("bray","unifrac","wunifrac","dpcoa","jsd"),"jaccard","philr")
  if (!distance_m %in% supported_methods){stop("Unsuported distance method \n\n")}
  if (distance_m %in% c("unifrac","wunifrac","dpcoa","jsd","bray")){
    d = phyloseq::distance(ps, distance_m)}
  if (distance_m %in% c("jaccard")){
    d = phyloseq::distance(ps, distance_m, binary=T)}
  if (distance_m=="philr"){
    suppressWarnings(
      ps.sub.philr <- philr(otu_table(ps) + 1, phy_tree(ps),part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt'))
    d = dist(ps.sub.philr, method="euclidean")}

  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  R2s <- pvals
  samples_used <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      samples_used[colnames(ps@sam_data)[i]] <- sum(!is.na(df[,i]))
      test.out <- adonis(as.matrix(d)[!is.na(df[,i]),!is.na(df[,i])] ~ get(colnames(ps@sam_data)[i])[!is.na(df[,i])], df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
      Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[1]
      R2s[colnames(ps@sam_data)[i]] <- test.out$aov.tab$R2[1]
    }, error=function(e){})
  }
  
  df.out <- cbind(samples_used[pvals!="NA"],pvals[pvals!="NA"],Fstat[pvals!="NA"],R2s[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"]), decreasing = T),]
  df.out <- setNames(data.frame(df.out), c("samples_used","p.value","F.model","R2"))
  df.out$F.model <- round(as.numeric(as.character(df.out$F.model)), digits = 4)
  df.out$p.value <- formatC(as.numeric(as.character(df.out$p.value)), format = "e", digits = 3)
  df.out$R2 <- round(as.numeric(as.character(df.out$R2)), digits = 4)

  return(df.out)
  }
