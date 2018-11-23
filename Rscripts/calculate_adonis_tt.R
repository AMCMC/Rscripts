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
calculate_adonis_bray_TT <- function(ps, report_top_taxa=5){
  require(phyloseq)
  require('vegan')
  #d = phyloseq::distance(ps, "bray")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  Fstat <- pvals
  R2s <- pvals
  tt <- pvals
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(ps@otu_table ~ get(colnames(ps@sam_data)[i]), method = "bray", df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
      Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[1]
      R2s[colnames(ps@sam_data)[i]] <- test.out$aov.tab$R2[1]
      if (report_top_taxa>0){
        tt[colnames(ps@sam_data)[i]] <- paste(names(sort(abs(test.out$coefficients[2,]),decreasing = T)[1:report_top_taxa]), collapse=";")
      }
    }, error=function(e){})
  }
  df <- cbind(pvals[pvals!="NA"],Fstat[pvals!="NA"],R2s[pvals!="NA"],tt[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"])),]
  setNames(data.frame(df), c("p.value","F.model","R2","toptaxa"))
}
plot_coef_top_adonis <- function(ps, var="Subject_ID", x=20, dist="bray"){
  require(vegan)
  require(ggplot2)
  df = as(sample_data(ps), "data.frame")
  perm_out <- adonis(ps@otu_table ~ get(var), method = dist, df)

  OTUS <- names(sort(abs(perm_out$coefficients[2,]), decreasing = T)[1:x])
  df <- data.frame(coef=perm_out$coefficients[2,][OTUS],
                   OTUS=OTUS,
                   tax_label=paste0(OTUS,";",apply(ps@tax_table[OTUS,2:7], 1, FUN = function(x) paste(x, collapse=";"))))
  #ggplot(df, aes(x=tax_label, y=coef, color=coef, fill=coef)) +
  #   geom_bar(stat="identity") +
  #   scale_fill_viridis_c() + 
  #   scale_colour_viridis_c() +
  #   coord_flip()
  df$tax_label <- factor(df$tax_label, df$tax_label[rev(order(df$coef))])
  
  #df <- df[rev(order(df$coef)),]
  p <- ggplot(df, aes(x=tax_label, y=coef, color=coef, fill=coef)) +   
    geom_bar(stat="identity") +
    scale_fill_viridis_c() + 
    scale_colour_viridis_c() +
    coord_flip()
  return(p)
}