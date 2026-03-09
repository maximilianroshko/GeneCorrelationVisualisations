library(ggplot2)
library(dplyr)
library(tidyr)
install.packages("ggrepel")
library(ggrepel)

rm(list = ls())

maindata <- read.csv("/Users/nj/Library/CloudStorage/OneDrive-QueenMary,UniversityofLondon/Year 2/Semester B/Rstudio project/B-cell lymphomas/Main dataset 1.csv")

# ── Step 1: Define cluster genes and stage order ──────────────────────────────
cluster_genes <- c("GAPDH", "DNMT1", "EZH2", "RAD51", "CCNB1", "POLH", "MKI67", "AURKA", "SMARCA4")  # replace with your cluster genes

stage_order <- c("GC LZ", "GC other", "GC DZ", "FL")  # <── using IndexedPheno values, GC other = IZ

# ── Step 2: Convert gene columns to numeric ───────────────────────────────────
cols <- which(colnames(maindata) == "AICDA"):which(colnames(maindata) == "TP53")
maindata[ , cols] <- lapply(maindata[ , cols], as.numeric)

# ── Step 3: Subset using IndexedPheno this time ───────────────────────────────
traj_data <- maindata %>%
  filter(IndexedPheno %in% stage_order,
         QCfilters == "pass") %>%
  select(IndexedPheno, all_of(cluster_genes)) %>%
  mutate(across(all_of(cluster_genes), as.numeric),
         IndexedPheno = factor(IndexedPheno, levels = stage_order))

# ── Quick sanity check — should show counts for all 4 stages ─────────────────
print(table(traj_data$IndexedPheno))

# ── Step 4: Mean expression per stage per gene ────────────────────────────────
traj_mean <- traj_data %>%
  group_by(IndexedPheno) %>%
  summarise(across(all_of(cluster_genes), mean, na.rm = TRUE)) %>%
  pivot_longer(-IndexedPheno, names_to = "Gene", values_to = "MeanExpression")

# ── Step 5: SD per stage per gene ─────────────────────────────────────────────
traj_sd <- traj_data %>%
  group_by(IndexedPheno) %>%
  summarise(across(all_of(cluster_genes), sd, na.rm = TRUE)) %>%
  pivot_longer(-IndexedPheno, names_to = "Gene", values_to = "SD")

# ── Step 6: Join and re-apply factor order ────────────────────────────────────
traj_summary <- left_join(traj_mean, traj_sd, by = c("IndexedPheno", "Gene")) %>%
  mutate(IndexedPheno = factor(IndexedPheno, levels = stage_order))

# ── Step 7: Labels at final stage (GC LZ) ────────────────────────────────────
label_data <- traj_summary %>% filter(IndexedPheno == "GC LZ")

# ── Step 8: Plot ──────────────────────────────────────────────────────────────
ggplot(traj_summary, aes(x = IndexedPheno, y = MeanExpression,
                         colour = Gene, group = Gene)) +
  
  geom_ribbon(aes(ymin = MeanExpression - SD,
                  ymax = MeanExpression + SD,
                  fill = Gene),
              alpha  = 0.12,
              colour = NA) +
  
  geom_line(linewidth = 1.3) +
  
  geom_point(size = 3.5, shape = 21,
             aes(fill = Gene), colour = "white", stroke = 0.6) +
  
  geom_text_repel(
    data          = label_data,
    aes(label     = Gene),
    nudge_x       = 0.25,
    direction     = "y",
    hjust         = 0,
    segment.color = "grey60",
    segment.size  = 0.4,
    size          = 3.8,
    fontface      = "bold"
  ) +
  
  geom_vline(xintercept = 1.5,
             linetype   = "dashed",
             colour     = "grey70",
             linewidth  = 0.5) +
  
  annotate("text", x = 1.52, y = Inf,
           label  = "GC entry →",
           hjust  = 0, vjust = 1.5,
           size   = 3, colour = "grey50") +
  
  # Rename the x axis tick for "GC other" so it reads more clearly on the plot
  scale_x_discrete(
    labels  = c("FL", "GC DZ", "GC IZ\n(other)", "GC LZ"),
    expand  = expansion(add = c(0.3, 0.8))
  ) +
  
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette   = "Set1") +
  
  labs(
    title    = "Gene Expression Across GC Trajectory",
    subtitle = "Shaded ribbon = ± 1 SD · FL used as pre-GC baseline · GC IZ = 'GC other' in IndexedPheno",
    x        = "Cell Stage",
    y        = "Mean Expression"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(colour = "grey50", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    legend.position  = "none",
    plot.margin      = margin(10, 40, 10, 10)
  )
