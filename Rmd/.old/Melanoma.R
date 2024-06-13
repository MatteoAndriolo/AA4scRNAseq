# %% Load necessary libraries
library(Seurat)
library(ggplot2)
library(here)
library(archetypes)
library(tidyr)
library(dplyr)
if (!require(kneedle)) {
  install.packages("quantmod")
}
devtools::install_github("etam4260/kneedle")
library(kneedle)

HVF <- TRUE
TEST <- TRUE

TEST_genes <- 300
TEST_samples <- 500

# %% Load Dataset
data_path <- "../data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt"
out_path <- "../out/Melanoma"

se <- read.table(data_path, header = TRUE)
se <- se[!duplicated(se[, 1]), ] # remove duplicated genes
rownames(se) <- se[, 1] # extract from matrix rownames
se <- se[, 2:ncol(se)] # elide rownames from gene expression matrix

metadata <- se[1:3, ] # extract metadata
metadata <- t(metadata) %>%
  data.frame() %>%
  mutate(across(where(is.character), as.numeric))

se <- se[4:nrow(se), ]
se <- se %>%
  data.frame() %>%
  mutate(across(where(is.character), as.numeric))
# -------- END BASIC SETUP
se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]

if (TEST) {
  tgenes <- min(TEST_genes, nrow(se))
  tsamples <- min(TEST_samples, ncol(se))
  metadata <- metadata[1:tsamples, ]
  se <- se[1:tgenes, 1:tsamples]
  rm(tgenes, tsamples)
  se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
}

se <- CreateSeuratObject(counts = se, meta.data = metadata)
se <- ScaleData(se, layer = "counts")
se <- FindVariableFeatures(se)

se <- RunPCA(se, features = VariableFeatures(se))
se <- RunUMAP(se, features = VariableFeatures(se))

# Visualize data
imgname <- sprintf("%s/PCA.png", out_path)
message(sprintf("Saving Imagine --- %s", imgname))
pcaplot <- PCAPlot(se)
# ggsave(imgname, plot=pcaplot)
print(pcaplot)

imgname <- sprintf("%s/UMAP.png", out_path)
message(sprintf("Saving Imagine --- %s", imgname))
umapplot <- UMAPPlot(se)
# ggsave(imgname, plot=umapplot)
print(umapplot)

imgname <- sprintf("%s/elbow.pdf", out_path)
message(sprintf("Saving Imagine --- %s", imgname))
elbowplot <- ElbowPlot(se)
# ggsave(imgname, plot=elbowplot)
print(elbowplot)

# Archetypes
k <- 10
k <- kneedle(elbowplot$data$dims, elbowplot$data$stdev)[1]

if (HVF) {
  m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
} else {
  m <- se@assays$RNA@layers$counts
}
m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
m <- as.matrix(m)

s <- timestamp()
a <- archetypes::archetypes(m, k = k, verbose = TRUE, maxIterations = 10)
e <- timestamp()

save(a, file = sprintf("%s/Archetypes_%2d.rds", out_path, k))

## Visualize Archetypes
# imgname <- sprintf("%s/Archetypes_%2d.png", out_path, k)
# plotarchetyps <- xyplot(a, m)
## ggsave(imgname, plot=plotarchetyps)
# print(plotarchetyps)

# Visualize UMAP archetypes
if (!require(umap)) {
  install.packages("umap")
}
library(umap)
library(cowplot)
umap_result <- umap(m)
umap_result <- umapplot[[1]]
umap_data <- as.data.frame(umap_result$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")
colnames(umap_data) <- c("umap_2", "umap_2")

# Get the archetype weights
weights <- coef(a)
weights <- as.data.frame(weights)

# Set a minimum threshold
threshold <- 0.2
weights[weights < threshold] <- 0

plot_list <- list()
# Plotting
for (i in 1:k) {
  umap_data$weight <- weights[, i]
  plot_title <- sprintf("UMAP Archetype %d", i)
  imgname <- sprintf("%s/UMAP_Archetype_%d.png", out_path, i)

  umap_plot <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = weight)) +
    geom_point(size = 1) +
    scale_color_gradient(low = "grey", high = "red") +
    ggtitle(plot_title) +
    labs(color = "Weight")
  umap_plot
  # ggsave(imgname, plot = umap_plot)

  # message(sprintf("Saving Image --- %s", imgname))
  plot_list[[i]] <- umap_plot
}
# Combine all plots into a single image
combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)

# Save the combined image
combined_imgname <- sprintf("%s/UMAP_Combined.png", out_path)
combined_plot
ggsave(combined_imgname, plot = combined_plot, width = 12, height = 8)
message(sprintf("Saving Combined Image --- %s", combined_imgname))



## Seurat
# DimPlot(se, reduction = "umap")
# DimPlot(se, reduction = "umap", group.by = "cell_type")
#
# se <- FindNeighbors(se, dims = 6:15)
# se <- FindClusters(se, resolution = 0.5)
