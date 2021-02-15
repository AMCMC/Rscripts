####  This script installs all required R packages from various sources ####
####  This script installs all required R packages from various sources ####

# cran is the main R package repository
# bioconda is another large repository with many usefull packages
# some people host small tools and scripts which can be usefull

#### cran packages ####

cran <- c("BiocManager",
          "reshape2",
          "ape",
          "devtools",
          "caret",
          "corrplot",
          "crayon",
          "dada2",
          "dplyr",
          "ggformula",
          "ggplot2",
          "ggpubr",
          "ggrepel",
          "ggsignif",
          "ggthemes",
          "ggbeeswarm",
          "ggtree",
          "gplots",
          "lme4",
          "lmerTest",
          "magick",
          "foreign",
          "kableExtra",
          "knitr",
          "matrixStats",
          "remotes",
          "permute",
          "phangorn",
          "philr",
          "phyloseq",
          "phytools",
          "picante",
          "rio",
          "xlsx",
          "vegan",
          "stringi",
          "stringr",
          "tibble",
          "tidyr",
          "tidyverse",
          "clusterSim",
          NULL)

install.packages(cran[!cran %in% installed.packages()])

#### bioc packages ####

bioc. <- c("phyloseq",
           "DESeq2",
           "dada2",
           "microbiome",
           "philr",
           "Biostrings",
           "DECIPHER",
           "decontam",
           "mixOmics",
           NULL)

BiocManager::install(bioc.[!bioc. %in% installed.packages()])

#### other sources ####

devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
remotes::install_github("jfq3/ggordiplots")
remotes::install_github("gavinsimpson/ggvegan")
devtools::install_github("gmteunisse/Fantaxtic")

# failed to install cran:
print(cran[!cran %in% installed.packages()])

# failed to install bioc
print(bioc.[!bioc. %in% installed.packages()])
