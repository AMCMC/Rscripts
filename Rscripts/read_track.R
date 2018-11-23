readtrack <- function(ps, rel=T){
  
  
    cbbPalette <- rev(c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999"))

track.df <- data.frame(row.names = rownames(ps@sam_data),
                       Sample = rownames(ps@sam_data),
                       Quality_filtered = ps@sam_data$readcount_input-ps@sam_data$readcount_filtered,
                       Merge_failed = ps@sam_data$readcount_filtered-ps@sam_data$readcount_merged,
                       Chimeras = ps@sam_data$readcount_merged-ps@sam_data$readcount_nonchim,
                       Mitochondria = ps@sam_data$readcount_mitochondria,
                       Chloroplast = ps@sam_data$readcount_chloroplast,
                       PCR_contaminations = ps@sam_data$readcount_PCR_contaminations,
                       Passed_processing = ps@sam_data$readcount_nonchim)

df.track.long <- melt(track.df,id.vars=c("Sample"))
df.track.long$read_count <- df.track.long$value

track.df.rel <- track.df
track.df.rel$Sample <- NULL
track.df.rel <- track.df.rel/rowSums(track.df.rel)
track.df.rel$Sample <- rownames(track.df.rel)
df.track.rel.long <- melt(track.df.rel,id.vars=c("Sample"))
df.track.rel.long$relative_contribution <- df.track.rel.long$value


if (rel!=T){
p <- ggplot(df.track.long,aes(Sample,read_count,fill=variable)) + geom_bar(stat="identity")
p <- p + theme_bw()
p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
p <- p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p <- p + guides(fill=guide_legend(ncol=1, reverse=F))
p <- p + scale_fill_manual(values=cbbPalette)
p <- p + labs(y="Total reads", title = "Overall sample read processing statistics")
} else{
p <- ggplot(df.track.rel.long,aes(Sample,relative_contribution,fill=variable)) + geom_bar(stat="identity")
p <- p + theme_bw()+ylab("Relative contribution")
p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
p <- p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p <- p + guides(fill=guide_legend(ncol=1, reverse=F))
p <- p + scale_fill_manual(values=cbbPalette)
p <- p + labs(y = "Relative contribution", title = "Overall sample read processing statistics")
}
return(p)

}
