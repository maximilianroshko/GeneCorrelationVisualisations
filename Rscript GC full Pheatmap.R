install.packages("tidyverse")
install.packages("remotes")
remotes::install_github("maximilianroshko/NMUltimate-R-Package")
library(ggplot2)
library(tidyverse)
install.packages("corrplot")
library(corrplot)

rm(list = ls())

maindata <- read.csv("C:/Users/bt24406/OneDrive - Queen Mary, University of London/Year 2/Semester B/Rstudio project/B-cell lymphomas/Main dataset 1.csv")

maindata <- read.csv("/Users/nj/Library/CloudStorage/OneDrive-QueenMary,UniversityofLondon/Year 2/Semester B/Rstudio project/B-cell lymphomas/Main dataset 1.csv")

cols <- which(colnames(maindata) == "AICDA"):which(colnames(maindata) == "TP53")
maindata[ , cols] <- lapply(maindata[ , cols], as.numeric)
gc_data <- subset(maindata, SimpleSortPheno == "GC" & QCfilters == "pass")

gene_cols <- which(colnames(gc_data) %in% c("AICDA", "TP53"))
gene_data <- gc_data[ , seq(min(gene_cols), max(gene_cols))]
gene_data <- as.data.frame(lapply(gene_data, as.numeric))

var_filter <- apply(gene_data, 2, var, na.rm = TRUE) > 0
var_filter[is.na(var_filter)] <- FALSE
gene_data <- gene_data[ , var_filter]

cor_full <- cor(gene_data, method = "spearman", use = "pairwise.complete.obs")
d_full <- as.dist(1 - cor_full)

pheatmap(
  cor_full,
  clustering_distance_rows = d_full,
  clustering_distance_cols = d_full,
  clustering_method = "average",
  color = colorRampPalette(c("blue", "white", "red"))(200),
  main = "All genes co-expression - GC DZ (Spearman)",
  fontsize = 6,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE)
