se <- as.matrix(obj@se@assays$RNA$counts)
aspe <- t(parameters(obj@archetypes$model))
rownames(aspe) <- rownames(se)
colnames(aspe) <- paste0("arch", 1:9)
newse <- cbind(se, aspe)
newse <- as(new, "dgCMatrix")

# Create a temporary Seurat object with the combined matrix
combined_obj <- Seurat::CreateSeuratObject(counts = newse)
combined_obj <- ScaleData(combined_obj, layer = "counts")
if (debug) message("DEBUG: obj_umapWithArchetypes | objdim ", dim(combined_obj)[[1]], " ", dim(combined_obj)[[2]])
combined_obj <- RunPCA(combined_obj, features = rownames(combined_obj))

# UMAP on combined matrix
combined_obj <- RunUMAP(combined_obj, dims = 1:20) # Adjust dimensions as needed
# Extract UMAP embeddings
umap_combined <- Embeddings(combined_obj, "umap")

# Separate the combined results for plotting
ct <- obj@se@meta.data$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK.
at <- colnames(aspe)
combined_umap_df <- data.frame(umap_combined)
# combined_umap_df$type <- rep(c("SE", "Archetype"), c(ncol(obj@se), ncol(aspe)))
# combined_umap_df$type <- rep(ct, c(ncol(obj@se), ncol(aspe)))
combined_umap_df$type <- c(ct, at)
shape_values <- setNames(c(16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26), c(unique(cell_types), archetype_labels))


# Plot for Approach 2
plot2 <- ggplot(combined_umap_df, aes(x = umap_1, y = umap_2, color = type)) +
  geom_point(data = subset(combined_umap_df, type == "SE"), aes(shape = type), alpha = 0.7, size = 1) +
  geom_point(data = subset(combined_umap_df, type == "Archetype"), aes(shape = type), size = 4) +
  scale_shape_manual(values = c(16, 17)) + # Use different shapes for SE and Archetypes
  theme_minimal() +
  labs(title = "UMAP Projection of Combined SE and Archetypes", x = "UMAP 1", y = "UMAP 2")

# plot2 <- ggplot(combined_umap_df, aes(x = umap_1, y = umap_2, color = type)) +
#   geom_point() +
#   theme_minimal() +
#   labs(title = "UMAP Projection of Combined SE and Archetypes", x = "UMAP 1", y = "UMAP 2")

# ggsave("umap_projection_2.png", plot2)
# if (debug) message("DEBUG: obj_umapArchetypes | Approach 2 plot saved")

obj@plots$umap_withArchetypes <- plot2

# obj <- obj_updateParams(obj, updateCurrent = TRUE, umap_threshold = treshold)


#   # Approach 2: Combine SE and Archetypes, Then Perform UMAP
# #se_matrix <- GetAssayData(obj@se, slot = "scale.data")
# #se_matrix <- obj@se

# # Transpose archetypes to match the SE samples
# archetypes_transposed <- t(parameters(obj@archetypes$model))
# if (debug) message("DEBUG: obj_umapArchetypes | archetypes_transposed dimension is ", dim(archetypes_transposed)[[1]], " ", dim(archetypes_transposed)[[2]])
# if (debug) message("DEBUG: obj_umapArchetypes | archetypes_transposed type is ", typeof(archetypes_transposed))

# # Combine SE and Archetypes
# #combined_matrix <- cbind(obj@se, archetypes_transposed)
# combined_matrix <- rbind(t(obj@se), parameters(obj@archetypes$model))

# # Create a temporary Seurat object with the combined matrix
# combined_obj <- CreateSeuratObject(counts = combined_matrix)
# combined_obj <- ScaleData(combined_obj)
# combined_obj <- RunPCA(combined_obj, features = rownames(combined_obj))

# # UMAP on combined matrix
# combined_obj <- RunUMAP(combined_obj, dims = 1:20)  # Adjust dimensions as needed

# # Extract UMAP embeddings
# umap_combined <- Embeddings(combined_obj, "umap")

# # Separate the combined results for plotting
# combined_umap_df <- data.frame(umap_combined)
# combined_umap_df$type <- rep(c('SE', 'Archetype'), c(ncol(obj@se), ncol(archetypes_transposed)))

# # Plot for Approach 2
# plot2 <- ggplot(combined_umap_df, aes(x = UMAP_1, y = UMAP_2, color = type)) +
#   geom_point() +
#   theme_minimal() +
#   labs(title = "UMAP Projection of Combined SE and Archetypes", x = "UMAP 1", y = "UMAP 2")

# #ggsave("umap_projection_2.png", plot2)
# # if (debug) message("DEBUG: obj_umapArchetypes | Approach 2 plot saved")

# obj@plots$umap_withArchetypes <- plot2

# obj <- obj_updateParams(obj, updateCurrent = TRUE, umap_threshold = treshold)

# # Approach 1: UMAP of SE and then Add Archetypes
# # Assuming obj@se has UMAP embeddings already computed
# umap_result <- UMAPPlot(obj@se, return.data = TRUE)
# umap_se <- as.data.frame(umap_result[, 1:2])
# colnames(umap_se) <- c("UMAP_1", "UMAP_2")

# # Transform archetypes using the same UMAP model
# umap_model <- obj@se@reductions$umap@misc$model
# umap_archetypes <- predict(umap_model, obj@archetypes$model$archetypes)

# # Combine results for plotting
# se_umap_df <- data.frame(umap_se, type = 'SE')
# archetypes_umap_df <- data.frame(umap_archetypes, type = 'Archetype')
# colnames(archetypes_umap_df) <- c("UMAP_1", "UMAP_2", "type")

# combined_df <- rbind(se_umap_df, archetypes_umap_df)

# # Plot for Approach 1
# plot1 <- ggplot(combined_df, aes(x = UMAP_1, y = UMAP_2, color = type)) +
#   geom_point() +
#   theme_minimal() +
#   labs(title = "UMAP Projection of SE with Archetypes", x = "UMAP 1", y = "UMAP 2")

# ggsave("umap_projection_1.png", plot1)
# if (debug) message("DEBUG: obj_umapArchetypes | Approach 1 plot saved")
return(obj)
