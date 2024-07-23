library(Rtsne)
library(plotly)
library(ggplot2)
library(cowplot)
library(networkD3)
if (!require(ggsankey)) devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)

set.seed(2024)
# Read files
se3D <- readRDS("data/Melanoma/GLYK.Rds")
se <- readRDS("data/Melanoma/GLYK.Rds")
se <- readRDS("data/Melanoma/se.Rds")
unique(se$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK.)

se$malignant.1.no.2.yes.0.unresolved. <- factor(
  se$malignant.1.no.2.yes.0.unresolved.,
  levels = c("1", "2", "0"),
  labels = c("non malignant", "malignant", "unresolved")
)
se$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK. <- factor(
  se$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK.,
  levels = c("0", "1", "2", "3", "4", "5", "6"),
  labels = c("?", "T", "B", "Macro", "Endo", "CAF", "NK")
)
table(se$malignant.1.no.2.yes.0.unresolved., se$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK.)
AA <- readRDS("data/Melanoma/GLYK_AA.Rds")

# Take one archetype
archetypes <- t(AA$aa.bests[["12"]]$BY)
num_archetypes <- ncol(archetypes)
common <- intersect(rownames(se), rownames(archetypes))
common

# Find closest points to archetypes
closest_points <- sapply(1:num_archetypes, function(i) {
  archetype <- archetypes[, i]
  distances <- apply(GetAssayData(se), 2, function(cell) sum((cell[common] - archetype)^2))
  which.min(distances)
})
closest_points
# Create a new ctype vector with archetypes marked
# new_ctype <- se$ctype
# levels(new_ctype)[closest_points] <- "Archetype"
# se$ctype <- factor(new_ctype, levels = c(levels(se$ctype), "Archetype"))

# Create a new ctype vector with archetypes marked
new_ctype <- as.character(se$ctype)
new_ctype[closest_points] <- "Archetype"
new_ctype <- factor(new_ctype, levels = c(levels(se$ctype), "Archetype"))

# Check the new ctype vector
table(new_ctype)


# Create Seurat object
newse <- CreateSeuratObject(counts = GetAssayData(se))
newse <- SetAssayData(newse, new.data = as.matrix(GetAssayData(se)), layer = "scale.data")
newse$ctype <- new_ctype

print(newse$ctype)
print(table(newse$ctype))

newse <- FindVariableFeatures(newse)
newse <- RunPCA(newse, features = VariableFeatures(newse))
newse <- RunUMAP(newse, features = VariableFeatures(newse))
DimPlot(newse, reduction = "umap")

# Run TSNE using Seurat's RunTSNE function
newse <- RunTSNE(newse, dims = 1:30, perplexity = 30, check_duplicates = FALSE, dim.embed = 3)
DimPlot(newse, reduction = "tsne")

tsne_3d <- as.data.frame(Embeddings(newse, reduction = "tsne"))
colnames(tsne_3d) <- c("X1", "X2", "X3")

tsne_3d$Label <- newse$ctype

# 2D Plots using ggplot2 with custom labels
labels_colors <- scales::hue_pal()(length(unique(tsne_3d$Label)))
archetype_color <- "black"

p1 <- ggplot(tsne_3d, aes(x = X1, y = X2, color = Label)) +
  geom_point(data = subset(tsne_3d, Label != "Archetype"), size = 1) +
  geom_point(data = subset(tsne_3d, Label == "Archetype"), color = "black", size = 4) +
  geom_text(
    data = subset(tsne_3d, Label == "Archetype"),
    aes(label = seq(1, length(closest_points))),
    color = "white", size = 3
  ) +
  scale_color_manual(values = c(labels_colors, archetype_color)) +
  theme_minimal() +
  ggtitle("X1 vs X2")

p2 <- ggplot(tsne_3d, aes(x = X1, y = X3, color = Label)) +
  geom_point(data = subset(tsne_3d, Label != "Archetype"), size = 1) +
  geom_point(data = subset(tsne_3d, Label == "Archetype"), color = "black", size = 4) +
  geom_text(
    data = subset(tsne_3d, Label == "Archetype"),
    aes(label = seq(1, length(closest_points))),
    color = "white", size = 3
  ) +
  scale_color_manual(values = c(labels_colors, archetype_color)) +
  theme_minimal() +
  ggtitle("X1 vs X3")

p3 <- ggplot(tsne_3d, aes(x = X2, y = X3, color = Label)) +
  geom_point(data = subset(tsne_3d, Label != "Archetype"), size = 1) +
  geom_point(data = subset(tsne_3d, Label == "Archetype"), color = "black", size = 4) +
  geom_text(
    data = subset(tsne_3d, Label == "Archetype"),
    aes(label = seq(1, length(closest_points))),
    color = "white", size = 3
  ) +
  scale_color_manual(values = c(labels_colors, archetype_color)) +
  theme_minimal() +
  ggtitle("X2 vs X3")

plot_3d <- plot_ly(tsne_3d, x = ~X1, y = ~X2, z = ~X3, color = ~Label, colors = c(labels_colors, archetype_color)) %>%
  add_markers(size = ~ ifelse(Label == "Archetype", 8, 3), text = ~ ifelse(Label == "Archetype", seq(1, length(closest_points)), "")) %>%
  layout(
    title = "3D t-SNE Visualization",
    scene = list(
      xaxis = list(title = "X1"),
      yaxis = list(title = "X2"),
      zaxis = list(title = "X3")
    )
  )

# Display the 3D plot
plot_3d

# Combine the plots into a single figure
combined_plot <- plot_grid(p1, p2, p3, ncol = 1)

# Save the combined plot as a PNG image
ggsave("combined_tsne_plot.png", combined_plot, width = 8, height = 12)

# Add labels for 3D plot
tsne_3d$size <- ifelse(tsne_3d$Label == "Archetype", 10, 5)
tsne_3d$color <- as.factor(ifelse(tsne_3d$Label == "Archetype", "black", tsne_3d$Label))
tsne_3d$GoldLabel <- sample(c(1, 4), num_points, replace = TRUE)

# 3D Plot using plotly
p3d <- plot_ly(tsne_3d,
  x = ~X1, y = ~X2, z = ~X3,
  type = "scatter3d", mode = "markers+text",
  marker = list(size = ~size, color = ~color, colorscale = "Viridis"),
  text = ifelse(tsne_3d$Label == "Archetype", rownames(tsne_3d), ""),
  textposition = "top center"
)

# SANKEY
# Assign ground truth labels randomly to the data
set.seed(123)
tsne_3d$GroundTruth <- sample(1:4, num_points, replace = TRUE)

# Assign numeric labels to factors with corresponding text labels
label_factors <- factor(tsne_3d$Label,
  levels = c(1:6, "Archetype"),
  labels = c("Type A", "Type B", "Type C", "Type D", "Type E", "Type F", "Archetype")
)

ground_truth_factors <- factor(tsne_3d$GroundTruth,
  levels = 1:4,
  labels = c("Class 1", "Class 2", "Class 3", "Class 4")
)

# Update the tsne_3d data frame with these factors
tsne_3d$LabelText <- label_factors
tsne_3d$GroundTruthText <- ground_truth_factors

# Create a data frame for the Sankey plot with text labels
sankey_data <- data.frame(
  From = tsne_3d$LabelText,
  To = tsne_3d$GroundTruthText
)

# Plot the Sankey plot with text labels
sankey_plot <- sankey_data %>%
  make_long(From, To) %>%
  ggplot(aes(
    x = x, next_x = next_x, node = node, next_node = next_node,
    fill = factor(node), label = node
  )) +
  geom_sankey(flow.alpha = 0.5, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white") +
  theme_sankey(base_size = 12) +
  labs(
    title = "Sankey Plot of Labels to Ground Truth",
    x = "Labels", y = "Ground Truth"
  )

# Save the Sankey plot as a PNG image
ggsave("sankey_plot.png", sankey_plot, width = 10, height = 6)
