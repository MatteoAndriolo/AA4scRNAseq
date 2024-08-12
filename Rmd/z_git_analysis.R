# z_git_analysis
library(Rtsne)
library(plotly)
library(ggplot2)
library(cowplot)
# library(networkD3)
if (!require(ggsankey)) devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(viridis)

set.seed(2024)

se <- readRDS("data/Melanoma/se.Rds")
name <- list(
  "GLYK",
  "MAPK",
  "CANCER",
  "MTOR",
  "TGF",
  "HVF"
)
pathw <- 1
se3D <- readRDS(paste0("data/Melanoma/", name[pathw], ".Rds"))
aa <- readRDS(paste0("data/Melanoma/", name[pathw], "_AA.Rds"))

addArchetypesToSeurat <- function(se, aa, k) {
  archetype <- t(aa$aa.bests[[k]]$BY)
  num_archetypes <- ncol(archetype)
  aa_clusters <- aa$aa.bests[[k]]$cluster.id
  aa_clusters
  aa_clusters_treshold <- aa$aa.bests[[k]]$cluster.id.treshold.5
  aa_clusters_treshold
  # weights dim = num_cells x num_archetypes
  weights <- aa$aa.bests[[k]]$A
  weights

  tempMat <- matrix(0, nrow = nrow(se), ncol = num_archetypes)
  rownames(tempMat) <- rownames(se)
  colnames(tempMat) <- paste("Archetype", 1:num_archetypes, sep = "_")

  common <- intersect(rownames(se), rownames(archetype))

  tempMat[common, ] <- archetype[common, ]

  tempMat <- as(tempMat, "dgCMatrix")

  tempMat <- cbind(GetAssayData(se), tempMat)

  newse <- CreateSeuratObject(counts = tempMat)
  newse <- SetAssayData(newse, new.data = as.matrix(tempMat), layer = "scale.data")
  newse <- FindVariableFeatures(newse)
  newse <- RunPCA(newse, features = VariableFeatures(newse))
  newse <- RunUMAP(newse, features = VariableFeatures(newse))
  newse <- RunTSNE(newse, dims = 1:30, perplexity = 30, check_duplicates = FALSE, dim.embed = 3)

  # Metadata
  tempType <- factor(
    c(as.character(se$ctype), rep("Archetype", num_archetypes)),
    levels = c(levels(se$ctype), "Archetype")
  )

  tempClusters <- factor(
    c(as.character(aa_clusters), rep("Archetype", num_archetypes)),
    labels = c(sort(unique(aa_clusters)), "Archetype")
  )
  tempClustersTreshold <- factor(
    c(as.character(aa_clusters_treshold), rep("Archetype", num_archetypes)),
    labels = c(sort(unique(aa_clusters_treshold)), "Archetype")
  )
  tempClustersTreshold

  tempWeights <- rbind(weights, diag(num_archetypes))

  newse$ctype <- tempType
  newse$ctype
  newse@misc$aa_clusters <- tempClusters
  newse@misc$aa_clusters
  newse@misc$aa_cluster_treshold.5 <- tempClustersTreshold
  newse@misc$aa_cluster_treshold.5
  # newse$aa_cluster_weight <- apply(tempWeights, 1,
  #     function(i){
  #         return(max(i))
  #     }
  # )
  newse@misc$weights <- tempWeights
  return(newse)
}

findClosestPoints <- function(se, aa, k) {
  archetype <- t(aa$aa.bests[[k]]$BY)
  aa_clusters <- aa$aa.bests[[k]]$cluster.id
  aa_clusters_treshold <- aa$aa.bests[[k]]$cluster.id.treshold.5
  num_archetypes <- ncol(archetype)
  weights <- as.data.frame(aa$aa.bests[[k]]$A)
  common <- intersect(rownames(se), rownames(archetype))

  small_se <- se[common, ]

  closestPoints <- sapply(1:num_archetypes, function(i) {
    distances <- apply(GetAssayData(small_se), 2, function(cell) sum((cell[common] - archetype[, i])^2))
    which.min(distances)
  })

  newCtype <- as.character(se$ctype)
  newCtype[closestPoints] <- "Archetype"
  newCtype <- factor(newCtype, levels = c(levels(se$ctype), "Archetype"))
  weights[closestPoints, ] <- diag(num_archetypes)
  rownames(weights) <- colnames(se)
  colnames(weights) <- paste("Archetype", 1:num_archetypes, sep = "_")

  se <- FindVariableFeatures(se)

  if (!"pca" %in% Reductions(se)) {
    se <- RunPCA(se, features = VariableFeatures(se))
  }
  se <- RunUMAP(se, features = VariableFeatures(se), dim.embed = 3)
  se <- RunTSNE(se, dims = 1:30, perplexity = 30, check_duplicates = FALSE, dim.embed = 3)

  se$ctype <- newCtype
  se@misc$aa_clusters <- aa_clusters
  se@misc$aa_cluster_treshold.5 <- aa$aa_bests[[k]]$cluster.id.treshold.5
  se@misc$weights <- weights

  return(se)
}


kappas <- list(
  "7",
  "12"
)
ik <- 2
k <- kappas[[ik]]
num_archetypes <- as.integer(kappas[[ik]])

################################################################################
# newse <- addArchetypesToSeurat(se, aa, "7")
newse <- findClosestPoints(se, aa, kappas[[ik]])
################################################################################

colors_types <- viridis(length(unique(newse$ctype)) - 1)
size_types <- 1
colors_Archetype <- "black"
size_Archetype <- 4
colors_text_Archetype <- c("white")
size_text_Archetype <- 3


plot_data <- as.data.frame(Embeddings(newse, reduction = "tsne"))
colnames(plot_data) <- c("X1", "X2", "X3")
plot_data$Label <- newse$ctype
plot_data$aa_clusters <- newse@misc$aa_clusters
plot_data$aa_clusters_treshold.5 <- newse@misc$aa_cluster_treshold.5

# 3D Plot
p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = Label)) +
  geom_point(data = subset(plot_data, Label != "Archetype"), size = size_types) +
  geom_point(data = subset(plot_data, Label == "Archetype"), color = "black", size = size_Archetype) +
  geom_text(
    data = subset(plot_data, Label == "Archetype"),
    aes(label = seq(1, num_archetypes)),
    color = colors_text_Archetype, size = size_text_Archetype
  ) +
  scale_color_manual(values = c(colors_types, colors_Archetype)) +
  theme_minimal()
p1

p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = Label)) +
  geom_point(data = subset(plot_data, Label != "Archetype"), size = size_types) +
  geom_point(data = subset(plot_data, Label == "Archetype"), color = "black", size = size_Archetype) +
  geom_text(
    data = subset(plot_data, Label == "Archetype"),
    aes(label = seq(1, num_archetypes)),
    color = colors_text_Archetype, size = size_text_Archetype
  ) +
  scale_color_manual(values = c(colors_types, colors_Archetype)) +
  theme_minimal()
p2

p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = Label)) +
  geom_point(data = subset(plot_data, Label != "Archetype"), size = size_types) +
  geom_point(data = subset(plot_data, Label == "Archetype"), color = "black", size = size_Archetype) +
  geom_text(
    data = subset(plot_data, Label == "Archetype"),
    aes(label = seq(1, num_archetypes)),
    color = colors_text_Archetype, size = size_text_Archetype
  ) +
  scale_color_manual(values = c(colors_types, colors_Archetype)) +
  theme_minimal()
p3
ggsave(
  file.path("data/Melanoma/", paste0(name[pathw], ".tsne.", ifelse(treshold > 0, paste(treshold, "."), ""), "X1.vs.X2.png")),
  p1,
  width = 8,
  height = 6
)
ggsave(
  file.path("data/Melanoma/", paste0(name[pathw], ".tsne.", ifelse(treshold > 0, paste(treshold, "."), ""), "X1.vs.X3.png")),
  p2,
  width = 8,
  height = 6
)
filename <- "data/Melanoma/test.png"
file.path("data/Melanoma/", paste0(name[pathw], ".tsne.", ifelse(treshold > 0, paste(treshold, "."), ""), "X2.vs.X3.png"))
ggsave(
  filename,
  plot = p3,
  width = 8,
  height = 6
)



plotIn3D <- function(plot_data, treshold = 0) {
}
treshold <- 0
plot_3d <- plotIn3D(plot_data, treshold)
plot_3d

# 2D plots of 3D

pairs <- list(
  c("X1", "X2"),
  c("X1", "X3"),
  c("X2", "X3")
)

treshold <- 0
for (p in pairs) {
  x <- p[1]
  y <- p[2]

  p2d <- ggplot(plot_data, aes(x = x, y = y, color = Label)) +
    geom_point(data = subset(plot_data, Label != "Archetype"), size = size_types) +
    geom_point(data = subset(plot_data, Label == "Archetype"), color = "black", size = size_Archetype) +
    geom_text(
      data = subset(plot_data, Label == "Archetype"),
      aes(label = seq(1, num_archetypes)),
      color = colors_text_Archetype, size = size_text_Archetype
    ) +
    scale_color_manual(values = c(colors_types, colors_Archetype)) +
    theme_minimal() +
    ggtitle(paste(name[pathw], " ", x, "vs", y, ifelse(treshold > 0, paste(" - treshold:", treshold), "")))

  print(p2d)
}

for (i in 1:length(pairs)) {
  t <- pairs[[i]]
}
combined_plot <- plot_grid(plotlist = plots, ncol = 1)
ggsave(
  paste0(name[pathw], "tsne.comb", ifelse(treshold > 0, paste(".", treshold), ""), ".png"),
)
