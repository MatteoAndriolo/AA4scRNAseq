library(Rtsne)
library(plotly)
library(ggplot2)
library(cowplot)
library(networkD3)
if (!require(ggsankey)) devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)

set.seed(2024)
# Read files
se <- readRDS("data/Melanoma/se.Rds")
{ # GLYc
  se3D <- readRDS("data/Melanoma/GLYK.Rds")
}
table(se$malignant.1.no.2.yes.0.unresolved., se$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK.)
AA <- readRDS("data/Melanoma/GLYK_AA.Rds")

# Take one archetype
archetypes <- t(AA$aa.bests[["7"]]$BY)

# Create matrix with archetypes
temp <- matrix(0, nrow = nrow(se), ncol = ncol(archetypes))
rownames(temp) <- rownames(se)
colnames(temp) <- paste("Archetype", 1:ncol(archetypes), sep = "_")
common <- intersect(rownames(se), rownames(archetypes))
temp[common, ] <- archetypes[common, ]
temp <- as(temp, "dgCMatrix")
temp <- cbind(GetAssayData(se), temp)

# Create Seurat object
newse <- CreateSeuratObject(counts = temp)
newse <- SetAssayData(newse, new.data = as.matrix(temp), layer = "scale.data")

arch_factor <- rep("Archetype", ncol(archetypes))
combined_ctype <- c(as.character(se$ctype), arch_factor)
newse$ctype <- factor(
  combined_ctype,
  levels = c(levels(se$ctype), "Archetype")
)
print(newse$ctype)
print(table(newse$ctype))

newse <- FindVariableFeatures(newse)
newse <- RunPCA(newse, features = VariableFeatures(newse))
newse <- RunUMAP(newse, features = VariableFeatures(newse))
DimPlot(newse, reduction = "umap")
# Run TSNE using Seurat's RunTSNE function
newse <- RunTSNE(newse, dims = 1:30, perplexity = 30, check_duplicates = FALSE, dim.embed = 3)
DimPlot(newse, reduction = "tsne")

# tsne_result <- Rtsne(t(GetAssayData(newse)), dims = 3, perplexity = 30, verbose = TRUE, max_iter = 500, partial_pca = TRUE)
# tsne_3d <- as.data.frame(tsne_result$Y)

# Perform t-SNE
# tsne_result <- Rtsne(data, dims = 3, perplexity = 30, verbose = TRUE, max_iter = 500)
tsne_3d <- as.data.frame(Embeddings(newse, reduction = "tsne"))
colnames(tsne_3d) <- c("X1", "X2", "X3")

# Assign random labels to the first n-4 points and "Archetype" to the last 4 points
num_points <- nrow(tsne_3d)
num_archetypes <- ncol(archetypes)
tsne_3d$Label <- newse$ctype
# plot <- ggplot(emb, aes(x = umap_1, y = umap_2, color = ctype)) +
#   geom_point(data = subset(emb, ctype != 99), size = 1) +
#   geom_point(data = subset(emb, ctype == 99), color = "black", size = 4) +
#   geom_text(
#     data = subset(emb, ctype == 99), aes(label = as.numeric(gsub("Archetype", "", rown))),
#     color = "white", size = 3
#   ) +
#   scale_color_manual(values = c(ctype_colors, rep(archetype_color, length(names_archetypes)))) +
#   theme_minimal() +
#   labs(title = "CellTypes and Archetypes found", x = "UMAP 1", y = "UMAP 2")
# 2D Plots using ggplot2 with custom labels
labels_colors <- scales::hue_pal()(length(unique(tsne_3d$Label)))
archetype_color <- "yellow"
p1 <- ggplot(tsne_3d, aes(x = X1, y = X2, color = Label)) +
  geom_point(data = subset(tsne_3d, Label != "Archetype"), size = 1) +
  geom_point(data = subset(tsne_3d, Label == "Archetype"), color = "black", size = 4) +
  # geom_point(aes(size = ifelse(Label == "Archetype", 5, 1)), show.legend = FALSE) +
  geom_text(
    data = subset(tsne_3d, Label == "Archetype"),
    aes(label = seq(1, num_archetypes)),
    color = "white", size = 3
  ) +
  scale_color_manual(values = c(labels_colors, rep(archetype_color, num_archetypes))) +
  theme_minimal() +
  ggtitle("X1 vs X2")

p2 <- ggplot(tsne_3d, aes(x = X1, y = X3, color = as.factor(Label))) +
  geom_point(aes(size = ifelse(Label == "Archetype", 5, 1)), show.legend = FALSE) +
  geom_text(
    data = subset(tsne_3d, Label == "Archetype"),
    aes(label = seq(num_points - 3, num_points)),
    color = "black", size = 5, vjust = -1
  ) +
  scale_color_manual(values = c(setNames(rep("grey", 6), 1:6), "black")) +
  theme_minimal() +
  ggtitle("X1 vs X3")

p3 <- ggplot(tsne_3d, aes(x = X2, y = X3, color = as.factor(Label))) +
  geom_point(aes(size = ifelse(Label == "Archetype", 5, 1)), show.legend = FALSE) +
  geom_text(
    data = subset(tsne_3d, Label == "Archetype"),
    aes(label = seq(num_points - 3, num_points)),
    color = "black", size = 5, vjust = -1
  ) +
  scale_color_manual(values = c(setNames(rep("grey", 6), 1:6), "black")) +
  theme_minimal() +
  ggtitle("X2 vs X3")

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

# Install and load the ggsankey package
install.packages("ggsankey")
library(ggsankey)

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
