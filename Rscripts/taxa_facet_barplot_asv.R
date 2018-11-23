taxa_facet_barplot_asv <-
  function(ps, Group, facet1, facet2, rank, lumpNA) {
    ps2 <- ps
    if (lumpNA) {
      if (dim(ps2@tax_table)[2] > 1) {
        if (sum(is.na(ps2@tax_table[, 2])) > 0) {
          ps2@tax_table[is.na(ps2@tax_table[, 2]), 2] <-
            paste(ps2@tax_table[is.na(ps2@tax_table[, 2]), 1], "_Unclassified", sep =
                    "")
        }
      }
      if (dim(ps2@tax_table)[2] > 2) {
        if (sum(is.na(ps2@tax_table[, 3])) > 0) {
          ps2@tax_table[is.na(ps2@tax_table[, 3]), 3] <-
            paste(ps2@tax_table[is.na(ps2@tax_table[, 3]), 2], "_Unclassified", sep =
                    "")
        }
      }
      if (dim(ps2@tax_table)[2] > 3) {
        if (sum(is.na(ps2@tax_table[, 4])) > 0) {
          ps2@tax_table[is.na(ps2@tax_table[, 4]), 4] <-
            paste(ps2@tax_table[is.na(ps2@tax_table[, 4]), 3], "_Unclassified", sep =
                    "")
        }
      }
      if (dim(ps2@tax_table)[2] > 4) {
        if (sum(is.na(ps2@tax_table[, 5])) > 0) {
          ps2@tax_table[is.na(ps2@tax_table[, 5]), 5] <-
            paste(ps2@tax_table[is.na(ps2@tax_table[, 5]), 4], "_Unclassified", sep =
                    "")
        }
      }
      if (dim(ps2@tax_table)[2] > 5) {
        if (sum(is.na(ps2@tax_table[, 6])) > 0) {
          ps2@tax_table[is.na(ps2@tax_table[, 6]), 6] <-
            paste(ps2@tax_table[is.na(ps2@tax_table[, 6]), 5], "_Unclassified", sep =
                    "")
        }
      }
      if (dim(ps2@tax_table)[2] > 6) {
        if (sum(is.na(ps2@tax_table[, 7])) > 0) {
          ps2@tax_table[is.na(ps2@tax_table[, 7]), 7] <-
            paste(ps2@tax_table[is.na(ps2@tax_table[, 7]), 6], "_Unclassified", sep =
                    "")
        }
      }
    }
    ps2@tax_table <-
      gsub("_Unclassified_.*Unclassified",
           "_Unclassified",
           ps2@tax_table)
    ps <- ps2
    colours <-
      c(
        "#F0A3FF",
        "#0075DC",
        "#993F00",
        "#4C005C",
        "#2BCE48",
        "#FFCC99",
        "#808080",
        "#94FFB5",
        "#8F7C00",
        "#9DCC00",
        "#C20088",
        "#003380",
        "#FFA405",
        "#FFA8BB",
        "#426600",
        "#FF0010",
        "#5EF1F2",
        "#00998F",
        "#740AFF",
        "#990000",
        "#FFFF00"
      )
    
    if (rank != "ASV") {
      ps <- tax_glom(ps, taxrank = rank, NArm = F)
    }
    
    if (taxa_are_rows(ps)) {
      abund_table <- t(ps@otu_table)
    } else {
      abund_table <- ps@otu_table
    }
    abund_table <- abund_table / rowSums(abund_table)
    abund_table <-
      abund_table[, order(colSums(abund_table), decreasing = TRUE)]
    N <- min(dim(abund_table)[2], length(colours) - 1)
    taxa_list <- rev(colnames(abund_table)[1:N])
    
    new_x <- NULL
    if (dim(abund_table)[2] > length(colours)) {
      new_x <-
        data.frame(Others = rowSums(abund_table[,!colnames(abund_table) %in% taxa_list]),
                   abund_table[, c(match(taxa_list, colnames(abund_table)))],
                   check.names = F)
    } else {
      new_x <- as.data.frame(abund_table[, taxa_list])
    }
    
    new_x$Sample <- rownames(new_x)
    if (exists("facet1")) {
      new_x$facet1 <-
        as.character(get(facet1, unclass(ps@sam_data[, facet1])))
    } else {
      print("Error: No Facet1 specified")
    }
    if (exists("facet2")) {
      new_x$facet2 <-
        as.character(get(facet2, unclass(ps@sam_data[, facet2])))
    } else {
      print("Error: No Facet2 specified")
    }
    if (exists("Group")) {
      new_x$Group <-
        as.factor(as.character(get(Group, unclass(ps@sam_data[, Group]))))
    } else {
      print("Error: No Group specified")
    }
    
    df <-
      reshape2::melt(new_x, value.name = "Value", variable.name = "Taxa")
    
    if (rank == "ASV") {
      #df$Taxa <- factor(df$Taxa,levels=rev(levels(df$Taxa)))
      for (Otu in grep("Other",
                       levels(df$Taxa),
                       invert = T,
                       value = T)) {
        levels(df$Taxa)[levels(df$Taxa) == Otu] <-
          paste(Otu, unname(ps@tax_table[Otu, 6]))
      }
    }
    if (rank != "ASV") {
      zotu = 22
      ps@tax_table[grep("Other",
                        levels(df$Taxa),
                        invert = T,
                        value = T), colnames(ps@tax_table) == rank]
      levels(df$Taxa)[match(grep(
        "Other",
        levels(df$Taxa),
        invert = T,
        value = T
      ), levels(df$Taxa))] <-
        ps@tax_table[grep("Other",
                          levels(df$Taxa),
                          invert = T,
                          value = T), colnames(ps@tax_table) == rank]
    }
    
    p <- ggplot(df, aes(Group, Value, fill = Taxa))
    p <- p + geom_bar(stat = "identity")
    p <-
      p + facet_grid(facet2 ~ facet1,
                     drop = TRUE,
                     scale = "free",
                     space = "free_x")
    if (dim(abund_table)[2] > length(colours)) {
      p <- p + scale_fill_manual(values = c(NA, rev(colours[0:(N)])))
    } else{
      p <- p + scale_fill_manual(values = c(rev(colours[0:(N)])))
    }
    p <- p + theme_bw() + ylab("Proportions")
    p <- p + scale_y_continuous(expand = c(0, 0)) +
      theme(strip.background = element_rect(fill = "gray85")) +
      theme(panel.spacing = unit(0.3, "lines"))
    p <-
      p + theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ))
    p <- p + guides(fill = guide_legend(ncol = 1, reverse = F))
    p <- p + theme(panel.spacing = unit(1, "lines"))
    p <-
      p + theme(text = element_text(size = 20),
                axis.text.x = element_text(angle = 90, hjust = 1))
    return(p)
  }