install.packages("igraph")
library(igraph)
install.packages("ggraph")
library(ggraph)
install.packages("threejs")
library(threejs)

rm(list = ls())

maindata <- read.csv("/Users/nj/Library/CloudStorage/OneDrive-QueenMary,UniversityofLondon/Year 2/Semester B/Rstudio project/B-cell lymphomas/Main dataset 1.csv")

cols <- which(colnames(maindata) == "AICDA"):which(colnames(maindata) == "TP53")
maindata[ , cols] <- lapply(maindata[ , cols], as.numeric)

gc_data <- subset(maindata, SimpleSortPheno == "GC" & QCfilters == "pass")

gene_cols <- which(colnames(gc_data) %in% c("AICDA", "TP53"))
gene_data <- gc_data[ , seq(min(gene_cols), max(gene_cols))]
gene_data <- gene_data[ , apply(gene_data, 2, var, na.rm = TRUE) > 0]

cor_full <- cor(gene_data, method = "spearman", use = "pairwise.complete.obs")

cor_pairs <- as.data.frame(as.table(cor_full))
colnames(cor_pairs) <- c("Gene1", "Gene2", "Correlation")
cor_pairs <- cor_pairs[cor_pairs$Gene1 != cor_pairs$Gene2, ]
cor_pairs <- cor_pairs[!duplicated(t(apply(cor_pairs[,1:2], 1, sort))), ]

# Asymmetric filter: strong positives and any meaningful negatives
edges <- cor_pairs[cor_pairs$Correlation >= 0.5 | cor_pairs$Correlation <= -0.2, ]
edges$Direction <- ifelse(edges$Correlation > 0, "Positive", "Negative")

cat("Positive edges (>= 0.5):", sum(edges$Direction == "Positive"), "\n")
cat("Negative edges (<= -0.2):", sum(edges$Direction == "Negative"), "\n")

all_genes <- colnames(gene_data)
gc_means  <- colMeans(gene_data, na.rm = TRUE)

g <- graph_from_data_frame(d = edges, directed = FALSE, vertices = data.frame(name = all_genes))
V(g)$mean_expression <- gc_means[V(g)$name]
g <- delete_vertices(g, V(g)[degree(g) == 0])

ggraph(g, layout = "fr") +
  geom_edge_link(aes(width  = abs(Correlation),
                     alpha  = abs(Correlation),
                     colour = Direction)) +
  scale_edge_colour_manual(values = c("Positive" = "#B2182B",
                                      "Negative" = "#2166AC"),
                           name   = "Correlation\nDirection") +
  scale_edge_width(range = c(0.3, 2.0), name = "Abs. Correlation\n(width)") +
  scale_edge_alpha(range = c(0.2, 0.9), name = "Abs. Correlation\n(alpha)") +
  geom_node_point(aes(size = mean_expression), colour = "#444444") +
  scale_size(range = c(2, 8), name = "Mean\nExpression") +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
  theme_void() +
  labs(title    = "Full Gene Co-expression Network (Spearman ≥ 0.5 | ≤ -0.2)",
       subtitle = "Red edges = positive correlation | Blue edges = negative correlation")

