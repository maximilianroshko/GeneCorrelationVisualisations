library(pheatmap)

rm(list = ls())

maindata <- read.csv("/Users/nj/Library/CloudStorage/OneDrive-QueenMary,UniversityofLondon/Year 2/Semester B/Rstudio project/B-cell lymphomas/Main dataset 1.csv")

cols <- which(colnames(maindata) == "AICDA"):which(colnames(maindata) == "TP53")
maindata[ , cols] <- lapply(maindata[ , cols], as.numeric)
gc_data <- subset(maindata, IndexedPheno == "GC LZ" & QCfilters == "pass")
gene_cols <- which(colnames(gc_data) %in% c("AICDA", "TP53"))
gene_data <- gc_data[ , seq(min(gene_cols), max(gene_cols))]
gene_data <- gene_data[ , apply(gene_data, 2, var, na.rm = TRUE) > 0]

cor_full <- cor(gene_data, method = "spearman", use = "pairwise.complete.obs")

mean_cor <- rowMeans(abs(cor_full))
top25_genes <- names(sort(mean_cor, decreasing = TRUE))[1:25]
cor_top25 <- cor_full[top25_genes, top25_genes]

d_top25 <- as.dist(1 - cor_top25)

pheatmap(
  cor_top25,
  clustering_distance_rows = d_top25,
  clustering_distance_cols = d_top25,
  clustering_method = "average",
  color = colorRampPalette(c("blue", "white", "red"))(200),
  main = "Top 25 co-expressed genes (Spearman)",
  fontsize = 10,
  border_color = "black")

