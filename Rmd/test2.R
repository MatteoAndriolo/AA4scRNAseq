source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")

obj <- readRDS("/app/out/Melanoma/Melanoma_files/latestUniqueH/MelanomaH_unique.rds")
# obj@archetypes$aa.kappas=list()
# saveRDS(obj, "/app/out/Melanoma/Melanoma_files/latestUniqueH/smallHunique.rds")
# obj <- readRDS("/app/out/Melanoma/Melanoma_files/latestUniqueH/smallHunique.rds")
classname <- "Melanoma"
class(obj) <- classname
# Extract and transpose archetype parameters
aspe <- t(parameters(obj@archetypes$model))
rownames(aspe) <- rownames(obj@se@assays$RNA$counts)
colnames(aspe) <- paste0("Archetype", 1:9)

# Combine the original matrix and archetypes
newse <- cbind(as.matrix(obj@se@assays$RNA$counts), aspe)
newse <- as(newse, "dgCMatrix")

# Create a temporary Seurat object with the combined matrix
combined_obj <- CreateSeuratObject(counts = newse)
combined_obj <- ScaleData(combined_obj, layer = "counts")
if (debug) message("DEBUG: obj_umapWithArchetypes | objdim ", dim(combined_obj)[[1]], " ", dim(combined_obj)[[2]])
combined_obj <- RunPCA(combined_obj, features = rownames(combined_obj))

# UMAP on combined matrix
combined_obj <- RunUMAP(combined_obj, dims = 1:20) # Adjust dimensions as needed

# Extract UMAP embeddings
umap_combined <- Embeddings(combined_obj, "umap")

# Separate the combined results for plotting
combined_umap_df <- data.frame(umap_combined)

# Add cell types for the original cells
cell_types <- obj@se@meta.data$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK.
archetype_labels <- colnames(aspe)
combined_umap_df$type <- c(cell_types, archetype_labels)

# Convert cell types to factors to ensure consistent ordering
combined_umap_df$type <- factor(combined_umap_df$type, levels = c(unique(cell_types), archetype_labels))

# Define color for archetypes and other points
archetype_color <- "yellow"
cell_type_colors <- scales::hue_pal()(length(unique(cell_types)))

# Plot UMAP with special points highlighted
plot2 <- ggplot(combined_umap_df, aes(x = umap_1, y = umap_2, color = type)) +
  geom_point(data = subset(combined_umap_df, !type %in% archetype_labels), size = 1) +
  geom_point(data = subset(combined_umap_df, type %in% archetype_labels), color = archetype_color, size = 4) +
  geom_text(
    data = subset(combined_umap_df, type %in% archetype_labels), aes(label = as.numeric(gsub("Archetype", "", type))),
    color = "black", size = 3
  ) + # vjust = -1.5, size = 3) +
  scale_color_manual(values = c(cell_type_colors, rep(archetype_color, length(archetype_labels)))) +
  theme_minimal() +
  labs(title = "UMAP Projection of Combined SE and Archetypes", x = "UMAP 1", y = "UMAP 2")

plot2

# Save or store the plot
# ggsave("umap_projection_2.png", plot2)
if (debug) message("DEBUG: obj_umapArchetypes | Approach 2 plot saved")

obj@plots$umap_withArchetypes <- plot2
