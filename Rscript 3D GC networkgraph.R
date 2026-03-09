install("igraph")
library(igraph)
install.packages("threejs")
library(threejs)
install.packages("htmlwidgets")
library(htmlwidgets)
install.packages("plotly")
library(plotly)

rm(list = ls())

#change this to read in the dataset.

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

edges <- cor_pairs[cor_pairs$Correlation >= 0.5 | cor_pairs$Correlation <= -0.2, ]
edges$Direction <- ifelse(edges$Correlation > 0, "Positive", "Negative")

cat("Positive edges (>= 0.5):", sum(edges$Direction == "Positive"), "\n")
cat("Negative edges (<= -0.2):", sum(edges$Direction == "Negative"), "\n")

all_genes <- colnames(gene_data)
gc_means  <- colMeans(gene_data, na.rm = TRUE)

g <- graph_from_data_frame(d = edges, directed = FALSE, vertices = data.frame(name = all_genes))
V(g)$mean_expression <- gc_means[V(g)$name]
g <- delete_vertices(g, V(g)[degree(g) == 0])

# 3D layout
set.seed(42)
layout_3d <- as.data.frame(layout_with_fr(g, dim = 3))
colnames(layout_3d) <- c("x", "y", "z")
layout_3d$name <- V(g)$name

# Node sizing
expr_vals  <- V(g)$mean_expression
node_sizes <- 4 + 12 * (expr_vals - min(expr_vals, na.rm = TRUE)) /
  (max(expr_vals, na.rm = TRUE) - min(expr_vals, na.rm = TRUE))

# Edge list
el        <- igraph::as_edgelist(g, names = TRUE)
edge_list <- data.frame(
  from        = el[, 1],
  to          = el[, 2],
  Correlation = E(g)$Correlation,
  Direction   = E(g)$Direction,
  stringsAsFactors = FALSE
)

# Edge trace builder
make_edge_traces <- function(edge_df, node_df, direction, colour) {
  sub <- edge_df[edge_df$Direction == direction, ]
  if (nrow(sub) == 0) return(NULL)
  
  abs_cor      <- abs(sub$Correlation)
  alpha_scaled <- 0.15 + 0.80 * (abs_cor - min(abs_cor)) /
    (max(abs_cor) - min(abs_cor) + 1e-9)
  width_scaled <- 0.5 + 3.5 * (abs_cor - min(abs_cor)) /
    (max(abs_cor) - min(abs_cor) + 1e-9)
  
  traces <- vector("list", nrow(sub))
  
  for (i in seq_len(nrow(sub))) {
    x0 <- node_df$x[node_df$name == sub$from[i]]
    y0 <- node_df$y[node_df$name == sub$from[i]]
    z0 <- node_df$z[node_df$name == sub$from[i]]
    x1 <- node_df$x[node_df$name == sub$to[i]]
    y1 <- node_df$y[node_df$name == sub$to[i]]
    z1 <- node_df$z[node_df$name == sub$to[i]]
    
    rgb_col  <- col2rgb(colour)
    rgba_col <- sprintf("rgba(%d,%d,%d,%.2f)",
                        rgb_col[1], rgb_col[2], rgb_col[3], alpha_scaled[i])
    
    traces[[i]] <- list(
      type        = "scatter3d",
      mode        = "lines",
      x           = c(x0, x1, NA),
      y           = c(y0, y1, NA),
      z           = c(z0, z1, NA),
      line        = list(color = rgba_col, width = width_scaled[i]),
      showlegend  = (i == 1),
      legendgroup = direction,
      name        = direction,
      hoverinfo   = "skip"
    )
  }
  return(traces)
}

pos_traces      <- make_edge_traces(edge_list, layout_3d, "Positive", "#B2182B")
neg_traces      <- make_edge_traces(edge_list, layout_3d, "Negative", "#2166AC")
all_edge_traces <- c(pos_traces, neg_traces)

# Node trace
node_trace <- list(
  type         = "scatter3d",
  mode         = "markers+text",
  x            = layout_3d$x,
  y            = layout_3d$y,
  z            = layout_3d$z,
  text         = layout_3d$name,
  textposition = "top center",
  textfont     = list(size = 9, color = "#222222"),
  hovertext    = paste0(layout_3d$name,
                        "<br>Mean Expression: ",
                        round(expr_vals, 3)),
  hoverinfo    = "text",
  marker       = list(
    size    = node_sizes,
    color   = "#444444",
    opacity = 0.85,
    line    = list(color = "white", width = 0.5)
  ),
  name        = "Nodes (size = Mean Expression)",
  showlegend  = TRUE,
  legendgroup = "nodes"
)

# Combine and plot
all_traces <- c(all_edge_traces, list(node_trace))

p <- plot_ly()
for (tr in all_traces) {
  p <- add_trace(p,
                 type         = tr$type,
                 mode         = tr$mode,
                 x            = tr$x,
                 y            = tr$y,
                 z            = tr$z,
                 line         = tr$line,
                 marker       = tr$marker,
                 text         = tr$text,
                 textposition = tr$textposition,
                 textfont     = tr$textfont,
                 hovertext    = tr$hovertext,
                 hoverinfo    = tr$hoverinfo,
                 name         = tr$name,
                 showlegend   = tr$showlegend,
                 legendgroup  = tr$legendgroup
  )
}

p <- layout(p,
            title = list(
              text = "3D Gene Co-expression Network (Spearman >= 0.5 | <= -0.2)",
              font = list(size = 14)
            ),
            scene = list(
              xaxis   = list(showticklabels = FALSE, title = "", showgrid = FALSE, zeroline = FALSE),
              yaxis   = list(showticklabels = FALSE, title = "", showgrid = FALSE, zeroline = FALSE),
              zaxis   = list(showticklabels = FALSE, title = "", showgrid = FALSE, zeroline = FALSE),
              bgcolor = "white"
            ),
            legend = list(
              title      = list(text = "<b>Legend</b>"),
              itemsizing = "constant",
              x          = 0.85,
              y          = 0.9
            ),
            paper_bgcolor = "white"
)

print(p)
