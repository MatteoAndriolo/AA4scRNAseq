# library(Rtsne)
# library(plotly)
# library(ggplot2)
# library(cowplot)
# if (!require(ggsankey)) devtools::install_github("davidsjoberg/ggsankey")
# if (!require(networkD3)) install.packages("networkD3")
# library(ggsankey)
# library(viridis)
# library(Seurat)
# library(scales)
# library(dplyr)
# if (!require(tidyverse)) install.packages("tidyverse")
source("Rmd/imports.R")
source("Rmd/classes.R")
source("Rmd/class_Melanoma.R")
source("Rmd/class_Mouse.R")
set.seed(2024)
NOT_FINAL <- TRUE
debug <- TRUE

obj <- new("Mouse")

if (NOT_FINAL) {
  plan("multicore", workers = 10)
  obj@params$hvf <- FALSE
  obj@params$test <- FALSE
  obj@params$pathw <- NULL
  temp_data_path <- file.path("data/MouseCortex/MouseCortex.RData")
}

obj <- obj_loadData(obj, data_path = temp_data_path)

if (NOT_FINAL) {
  obj@params$hvf <- FALSE
  obj@params$pathw <- 1

  obj@params$out_path <- "outPar/Mouse/0813_0758/FS1_1738962/"
  obj@params$path_figures <- file.path(obj@params$out_path, "figures")
  obj@params$path_outdata <- file.path(obj@params$out_path, "data")
  insefile <- file.path(obj@params$path_outdata, "1.Rds")
  inaafile <- file.path(obj@params$path_outdata, "archetypes.Rds")
  inmdfile <- file.path(obj@params$path_outdata, "metadata.Rds")


  obj@se <- readRDS(insefile)
  obj@archetypes <- readRDS(inaafile)
  obj@other <- readRDS(inmdfile)
}

##################################################
if (obj@params$hvf) {
  obj@other$namePathw <- "HVF"
} else {
  obj@other$namePathw <- list(
    "GLYK",
    "MAPK",
    "CANCER",
    "MTOR",
    "TGF"
  )[[obj@params$pathw]]
}
obj@other$treshold <- 0.5



# Function Definitions ---------------------------------------------------------
addArchetypesToSeurat <- function(se, aa, k) {
  archetype <- t(aa$aa.bests[[k]]$BY)
  num_archetypes <- ncol(archetype)
  aa_clusters <- aa$aa.bests[[k]]$cluster.id
  # aa_clusters
  # TODO remove comment from treshold once generated files in main
  # aa_clusters_treshold <- aa$aa.bests[[k]]$cluster.id.treshold
  # aa_clusters_treshold

  weights <- aa$aa.bests[[k]]$A
  weights

  tempMat <- matrix(0, nrow = nrow(se), ncol = num_archetypes)
  rownames(tempMat) <- rownames(se)
  colnames(tempMat) <- paste0("Archetype", 1:num_archetypes, sep = "_")

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

# FUnction definition ---------------------------------------------------------
findClosestPoints <- function(se, se3D, aa, k) {
  if (length(colnames(se)) != length(colnames(se3D))) {
    se <- se[, colnames(se3D)]
  }
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

  while (any(duplicated(closestPoints))) {
    firstdup <- which(duplicated(closestPoints))[1]
    firstduppoints <- closestPoints[5]
    distances <- apply(GetAssayData(small_se), 2, function(cell) sum((cell[common] - archetype[, 5])^2))
    while (which.min(distances) %in% closestPoints) {
      distances[which.min(distances)] <- Inf
    }
    closestPoints[firstdup] <- which.min(distances)
  }

  message("Closest points are: ", print(closestPoints))

  newCtype <- as.character(se$ctype)
  newCtype[closestPoints] <- "Archetype"
  newCtype <- factor(newCtype, levels = c(levels(se$ctype), "Archetype"))
  weights[closestPoints, ] <- diag(num_archetypes)
  rownames(weights) <- colnames(se)
  colnames(weights) <- paste0("Archetype", 1:num_archetypes, sep = "_")

  se <- FindVariableFeatures(se)

  if (!"pca" %in% Reductions(se)) {
    se <- RunPCA(se, features = VariableFeatures(se))
  }
  se <- RunUMAP(se, features = VariableFeatures(se), dim.embed = 3, n.components = 3)
  se <- RunTSNE(se, dims = 1:30, perplexity = 30, check_duplicates = FALSE, dim.embed = 3)

  se$ctype <- newCtype
  se@misc$aa_clusters <- aa_clusters
  se@misc$aa_cluster_treshold.5 <- aa$aa.bests[[k]]$cluster.id.treshold.5
  se@misc$weights <- weights

  return(se)
}

# Main -------------------------------------------------------------------------
if (FALSE) {
  k <- "7"
  i <- 1
  red <- "tsne"
}

for (k in c("7", "12")) {
  message("Using k=", k)
  num_archetypes <- as.integer(k)

  message("Pathw ", namePathw, " with k ", k)
  if (!dir.exists(obj@params$out_path)) dir.create(obj@params$out_path)

  # se3D <- readRDS(file.path(in_path, paste0(namePathw, ".Rds")))
  # aa <- readRDS(file.path(in_path, paste0(namePathw, "_AA.Rds")))

  ################################################################################
  # newse <- addArchetypesToSeurat(se, aa, "7")
  newse <- findClosestPoints(obj@se.org, obj@se, obj@archetypes, k)
  ################################################################################
  newse$ctype <- factor(newse$ctype)
  newse$aaclusters <- factor(obj@archetypes$aa.bests[[k]]$cluster.id)
  t <- obj@archetypes$aa.bests[[k]]$cluster.id.treshold
  t[t == "NA"] <- "NotAssigned"
  newse$aaclusters.treshold <- factor(t)
  newse$aaweights <- apply(obj@archetypes$aa.bests[["7"]]$A, 1, max)



  if (NOT_FINAL) {
    t <- obj@archetypes$aa.bests[[k]]$cluster.id
    t[t == "NA"] <- "1"
    newse$aaclusters <- factor(t)
  }
  if (NOT_FINAL) {
    t <- obj@archetypes$aa.bests[[k]]$cluster.id
    t[t == "NA"] <- "NotAssigned"
    newse$aaclusters.treshold <- factor(t)
  }
  all_colors <- rainbow(
    length(unique(newse$ctype)) + num_archetypes + 1
  )

  levels(newse$aaclusters) <- c(levels(newse$aaclusters), "Archetype")
  newse$aaclusters[newse$ctype == "Archetype"] <- "Archetype"
  newse$aaclusters <- factor(newse$aaclusters)

  levels(newse$aaclusters.treshold) <- c(levels(newse$aaclusters.treshold), "Archetype")
  newse$aaclusters.treshold[newse$ctype == "Archetype"] <- "Archetype"
  newse$aaclusters.treshold <- factor(newse$aaclusters.treshold)

  # CTYPES styling
  all_colorsCTypes <- all_colors[1:(length(unique(newse$ctype)) - 1)]
  colorMapCTypes <- setNames(
    c(all_colorsCTypes, "black"),
    levels(newse$ctype)
  )

  sizeMapCTypes <- setNames(
    c(rep(1, length(levels(newse$ctype)) - 1), 4),
    levels(newse$ctype)
  )
  size_text_Archetype <- 3



  # ARCHETYPES colors
  all_colorsArchetypes <- all_colors[
    (length(unique(newse$ctype))):(length(unique(newse$ctype)) + num_archetypes - 1)
  ]
  colorMapArchetypes <- setNames(
    c(all_colorsArchetypes, "grey", "black"),
    c(paste0("Archetype", 1:num_archetypes, sep = "_"), "NotAssigned", "Archetype")
  )
  colorMapArchetypes <- setNames(
    c(all_colorsArchetypes, "grey", "black"),
    c(as.character(1:num_archetypes), "NotAssigned", "Archetype")
  )
  sizeMapArchetypes <- setNames(
    c(rep(1, num_archetypes + 1), 4),
    c(as.character(1:num_archetypes), "NotAssigned", "Archetype")
  )

  colors_text_Archetype <- c("white")

  # Other Styling
  if (class(obj) == "Mouse") {
    newse$Time_points <- factor(newse$Time_points)
    shapeMapTimePoints <- setNames(
      c(15, 16, 17, 18),
      levels(newse$Time_points)
    )
  }
  if (class(obj) == "Melanoma") {
    newse$malignant <- factor(newse$malignant)
    shapeMapMalignant <- setNames(
      c(16, 17, 18),
      levels(newse$malignant)
    )
  }


  # for (red in c("umap", "tsne")) {
  plot_data <- as.data.frame(Embeddings(newse, reduction = red))
  colnames(plot_data) <- c("X1", "X2", "X3")
  plot_data <- cbind(plot_data, as.data.frame(newse@meta.data))
  for (red in c("tsne")) {
    message("reduction ", red)

    ##################################################
    # CELL TYPES
    ##################################################
    # 3D Plot
    p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = ctype, size = ctype)) +
      geom_point(data = subset(plot_data, ctype != "Archetype")) +
      geom_point(data = subset(plot_data, ctype == "Archetype")) +
      geom_text(
        data = subset(plot_data, ctype == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapCTypes) +
      scale_size_manual(values = sizeMapCTypes) +
      labs(color = "Cell types", size = "Cell types") +
      theme_classic()

    # p1

    p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = ctype, size = ctype)) +
      geom_point(data = subset(plot_data, ctype != "Archetype")) +
      geom_point(data = subset(plot_data, ctype == "Archetype")) +
      geom_text(
        data = subset(plot_data, ctype == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapCTypes) +
      scale_size_manual(values = sizeMapCTypes) +
      labs(color = "Cell types", size = "Cell types", ) +
      theme_classic()
    # p2

    p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = ctype, size = ctype)) +
      geom_point(aes(color = ctype, size = ctype),
        data = subset(plot_data, ctype != "Archetype")
      ) +
      geom_point(aes(color = ctype, size = ctype),
        data = subset(plot_data, ctype == "Archetype")
      ) +
      geom_text(
        data = subset(plot_data, ctype == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapCTypes) +
      scale_size_manual(values = sizeMapCTypes) +
      labs(color = "Cell types", size = "Cell types") +
      theme_classic()
    # p3

    combined_plot <- plot_grid(p1, p2, p3, ncol = 1)
    combined_plot

    prefixName <- paste(namePathw, k, red, sep = ".")
    ggsave(
      file.path(obj@params$out_path, paste0(prefixName, ".ct.X1.vs.X2.png")),
      p1,
      width = 8,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste0(prefixName, ".ct.X1.vs.X3.png")),
      p2,
      width = 8,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste0(prefixName, ".ct.X2.vs.X3.png")),
      plot = p3,
      width = 8,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste0(prefixName, ".ct.png")),
      combined_plot,
      width = 8,
      height = 18
    )


    if (class(obj) == "Mouse") { # TIME
      p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = ctype, size = ctype, shape = Time_points)) +
        geom_point(data = subset(plot_data, ctype != "Archetype")) +
        geom_point(data = subset(plot_data, ctype == "Archetype")) +
        geom_text(
          data = subset(plot_data, ctype == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapCTypes) +
        scale_size_manual(values = sizeMapCTypes) +
        labs(color = "Cell types", size = "Cell types") +
        theme_classic()

      # p1

      p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = ctype, size = ctype, shape = Time_points)) +
        geom_point(data = subset(plot_data, ctype != "Archetype")) +
        geom_point(data = subset(plot_data, ctype == "Archetype")) +
        geom_text(
          data = subset(plot_data, ctype == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapCTypes) +
        scale_size_manual(values = sizeMapCTypes) +
        labs(color = "Cell types", size = "Cell types", ) +
        theme_classic()
      # p2

      p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = ctype, size = ctype, shape = Time_points)) +
        geom_point(aes(color = ctype, size = ctype),
          data = subset(plot_data, ctype != "Archetype")
        ) +
        geom_point(aes(color = ctype, size = ctype),
          data = subset(plot_data, ctype == "Archetype")
        ) +
        geom_text(
          data = subset(plot_data, ctype == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapCTypes) +
        scale_size_manual(values = sizeMapCTypes) +
        labs(color = "Cell types", size = "Cell types") +
        theme_classic()
      # p3

      combined_plot <- plot_grid(p1, p2, p3, ncol = 1)
      combined_plot

      prefixName <- paste(namePathw, k, red, "time", sep = ".")
      ggsave(
        file.path(obj@params$out_path, paste0(prefixName, ".ct.X1.vs.X2.png")),
        p1,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste0(prefixName, ".ct.X1.vs.X3.png")),
        p2,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste0(prefixName, ".ct.X2.vs.X3.png")),
        plot = p3,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste0(prefixName, ".ct.png")),
        combined_plot,
        width = 8,
        height = 18
      )
    }


    ##################################################
    # Archetypes colored
    ##################################################
    # No Treshold
    p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = aaclusters, size = aaclusters)) +
      geom_point(data = subset(plot_data, aaclusters != "Archetype")) +
      geom_point(data = subset(plot_data, aaclusters == "Archetype")) +
      geom_text(
        data = subset(plot_data, aaclusters == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapArchetypes) +
      scale_size_manual(values = sizeMapArchetypes) +
      labs(color = "", size = "") +
      theme_classic()

    # p1

    p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = aaclusters, size = aaclusters)) +
      geom_point(aes(color = aaclusters, size = aaclusters),
        data = subset(plot_data, aaclusters != "Archetype")
      ) +
      geom_point(aes(color = aaclusters, size = aaclusters),
        data = subset(plot_data, aaclusters == "Archetype")
      ) +
      geom_text(
        data = subset(plot_data, aaclusters == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapArchetypes) +
      scale_size_manual(values = sizeMapArchetypes) +
      labs(color = "", size = "") +
      theme_classic()
    # p2

    p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = aaclusters, size = aaclusters)) +
      geom_point(data = subset(plot_data, aaclusters != "Archetype")) +
      geom_point(data = subset(plot_data, aaclusters == "Archetype")) +
      geom_text(
        data = subset(plot_data, aaclusters == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapArchetypes) +
      scale_size_manual(values = sizeMapArchetypes) +
      labs(color = "", size = "") +
      theme_classic()
    # p3

    combined_plot <- plot_grid(p1, p2, p3, ncol = 1)
    combined_plot

    prefixName <- paste(namePathw, k, red, ifelse(obj@other$treshold > 0, "th", ""), sep = ".")
    ggsave(
      file.path(obj@params$out_path, paste(prefixName, ".aa.X1.vs.X2.png")),
      p1,
      width = 9,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste(prefixName, ".aa.X1.vs.X3.png")),
      p2,
      width = 9,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste(prefixName, ".aa.X2.vs.X3.png")),
      plot = p3,
      width = 9,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste0(prefixName, ".aa.comb.png")),
      combined_plot,
      width = 9,
      height = 18
    )
    # Yes Treshold
    p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = aaclusters.treshold, size = aaclusters.treshold)) +
      geom_point(data = subset(plot_data, aaclusters.treshold != "Archetype")) +
      geom_point(data = subset(plot_data, aaclusters.treshold == "Archetype")) +
      geom_text(
        data = subset(plot_data, aaclusters.treshold == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapArchetypes) +
      scale_size_manual(values = sizeMapArchetypes) +
      theme_classic()

    # p1

    p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = aaclusters.treshold, size = aaclusters.treshold)) +
      geom_point(data = subset(plot_data, aaclusters.treshold != "Archetype")) +
      geom_point(data = subset(plot_data, aaclusters.treshold == "Archetype")) +
      geom_text(
        data = subset(plot_data, aaclusters.treshold == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapArchetypes) +
      scale_size_manual(values = sizeMapArchetypes) +
      theme_classic()
    # p2

    p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = aaclusters.treshold, size = aaclusters.treshold)) +
      geom_point(data = subset(plot_data, aaclusters.treshold != "Archetype")) +
      geom_point(data = subset(plot_data, aaclusters.treshold == "Archetype")) +
      geom_text(
        data = subset(plot_data, aaclusters.treshold == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = colorMapArchetypes) +
      scale_size_manual(values = sizeMapArchetypes) +
      theme_classic()
    # p3

    combined_plot <- plot_grid(p1, p2, p3, ncol = 1)

    prefixName <- paste(namePathw, k, red, ifelse(obj@other$treshold > 0, "th", ""), sep = ".")
    ggsave(
      file.path(obj@params$out_path, paste(prefixName, ".aa.X1.vs.X2.png")),
      p1,
      width = 8,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste0(prefixName, ".aa.X1.vs.X3.png")),
      p2,
      width = 8,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste(prefixName, ".aa.X2.vs.X3.png")),
      plot = p3,
      width = 8,
      height = 6
    )
    ggsave(
      file.path(obj@params$out_path, paste0(prefixName, ".aa.comb.png")),
      combined_plot,
      width = 8,
      height = 18
    )

    if (class(obj) == "Mouse") {
      # No Treshold
      p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = aaclusters, size = aaclusters, shape=Time_points)) +
        geom_point(data = subset(plot_data, aaclusters != "Archetype")) +
        geom_point(data = subset(plot_data, aaclusters == "Archetype")) +
        geom_text(
          data = subset(plot_data, aaclusters == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapArchetypes) +
        scale_size_manual(values = sizeMapArchetypes) +
        labs(color = "", size = "") +
        theme_classic()

      # p1

      p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = aaclusters, size = aaclusters, shape=Time_points)) +
        geom_point(aes(color = aaclusters, size = aaclusters),
          data = subset(plot_data, aaclusters != "Archetype")
        ) +
        geom_point(aes(color = aaclusters, size = aaclusters),
          data = subset(plot_data, aaclusters == "Archetype")
        ) +
        geom_text(
          data = subset(plot_data, aaclusters == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapArchetypes) +
        scale_size_manual(values = sizeMapArchetypes) +
        labs(color = "", size = "") +
        theme_classic()
      # p2

      p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = aaclusters, size = aaclusters, shape=Time_points)) +
        geom_point(data = subset(plot_data, aaclusters != "Archetype")) +
        geom_point(data = subset(plot_data, aaclusters == "Archetype")) +
        geom_text(
          data = subset(plot_data, aaclusters == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapArchetypes) +
        scale_size_manual(values = sizeMapArchetypes) +
        labs(color = "", size = "") +
        theme_classic()
      # p3

      combined_plot <- plot_grid(p1, p2, p3, ncol = 1)
      combined_plot

      prefixName <- paste(namePathw, k, red,"time", ifelse(obj@other$treshold > 0, "th", ""), sep = ".")
      ggsave(
        file.path(obj@params$out_path, paste(prefixName, ".aa.X1.vs.X2.png")),
        p1,
        width = 9,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste(prefixName, ".aa.X1.vs.X3.png")),
        p2,
        width = 9,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste(prefixName, ".aa.X2.vs.X3.png")),
        plot = p3,
        width = 9,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste0(prefixName, ".aa.comb.png")),
        combined_plot,
        width = 9,
        height = 18
      )
      # Yes Treshold
      p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = aaclusters.treshold, size = aaclusters.treshold, shape=Time_points)) +
        geom_point(data = subset(plot_data, aaclusters.treshold != "Archetype")) +
        geom_point(data = subset(plot_data, aaclusters.treshold == "Archetype")) +
        geom_text(
          data = subset(plot_data, aaclusters.treshold == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapArchetypes) +
        scale_size_manual(values = sizeMapArchetypes) +
        theme_classic()

      # p1

      p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = aaclusters.treshold, size = aaclusters.treshold, shape=Time_points)) +
        geom_point(data = subset(plot_data, aaclusters.treshold != "Archetype")) +
        geom_point(data = subset(plot_data, aaclusters.treshold == "Archetype")) +
        geom_text(
          data = subset(plot_data, aaclusters.treshold == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapArchetypes) +
        scale_size_manual(values = sizeMapArchetypes) +
        theme_classic()
      # p2

      p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = aaclusters.treshold, size = aaclusters.treshold, shape=Time_points)) +
        geom_point(data = subset(plot_data, aaclusters.treshold != "Archetype")) +
        geom_point(data = subset(plot_data, aaclusters.treshold == "Archetype")) +
        geom_text(
          data = subset(plot_data, aaclusters.treshold == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = colorMapArchetypes) +
        scale_size_manual(values = sizeMapArchetypes) +
        theme_classic()
      # p3

      combined_plot <- plot_grid(p1, p2, p3, ncol = 1)

      prefixName <- paste(namePathw, k, red,"time", ifelse(obj@other$treshold > 0, "th", ""), sep = ".")
      ggsave(
        file.path(obj@params$out_path, paste(prefixName, ".aa.X1.vs.X2.png")),
        p1,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste0(prefixName, ".aa.X1.vs.X3.png")),
        p2,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste(prefixName, ".aa.X2.vs.X3.png")),
        plot = p3,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(obj@params$out_path, paste0(prefixName, ".aa.comb.png")),
        combined_plot,
        width = 8,
        height = 18
      )
    }
  }
  # end red
  ggplot(plot_data, aes(x = aaweights)) +
    geom_histogram(
      binwidth = 0.05, # Adjust binwidth to change the number of bins
      fill = "blue",
      color = "black"
    ) +
    labs(
      x = "Weights",
      y = "Frequency"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    theme_classic()
  ggsave(
    file.path(obj@params$out_path, paste0(namePathw, k, "weights.png")),
    width = 8,
    height = 6
  )

  # 3D Plot
  p3d <- plot_ly(
    plot_data,
    x = ~X1, y = ~X2, z = ~X3,
    color = ~aaclusters.treshold,
    colors = colorMapArchetypes,
    size = sizeMapArchetypes[plot_data$aaclusters.treshold],
    type = "scatter3d",
    text = ~ ifelse(aaclusters == "Archetype", seq(1, num_archetypes), ""),
    textposition = "top center",
    mode = "markers",
    marker = list(
      opacity = 0.9
    )
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "X1"),
        yaxis = list(title = "X2"),
        zaxis = list(title = "X3")
      )
    )
  # p3d

  # SANKEY PLOT
  # Assuming newse$ctype and newse@misc$aa_clusters are vectors of the same length
  for (treshold in c(TRUE, FALSE)) {
    # taa_clusters <- as.character(newse$aaclusters)
    # taa_clusters <- plot_data$aaclusters
    # if (treshold) {
    #   taa_clusters <- plot_data$aaclusters.treshold
    # }

    # Create a data frame with the required columns
    if (treshold) {
      data <- as.data.frame(list(
        type = plot_data$ctype,
        archetype = plot_data$aaclusters.treshold
      ))
    } else {
      data <- as.data.frame(list(
        type = plot_data$ctype,
        archetype = plot_data$aaclusters
      ))
    }

    data$archetype <- factor(
      data$archetype,
      levels = c("1", "2", "3", "4", "5", "6", "7", "NotAssigned", "Archetype"),
      labels = c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "NotAssigned", "Archetype")
    )

    d <- data %>%
      make_long(colnames(data)) %>%
      filter(node != "Archetype") # %>% filter(next_node != "Archetype")
    d


    colorMapArchetypesSankey <- setNames(
      colorMapArchetypes,
      c(paste0("A", 1:num_archetypes, sep = ""), "NotAssigned", "Archetype")
    )

    pl <- ggplot(d, aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = factor(node),
      label = node
    )) +
      geom_sankey(
        flow.alpha = 0.5,
        node.color = "black",
        show.legend = FALSE
      ) +
      geom_sankey_label(size = 3, color = "black", fill = "white") +
      scale_fill_manual(
        values = c(colorMapCTypes, colorMapArchetypesSankey)
      ) +
      scale_x_discrete(
        labels = c("type" = "Cell Type", "archetype" = "Archetype")
      ) +
      theme_alluvial() +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    pl

    prefixName <- paste(namePathw, k, ifelse(treshold > 0, "th", ""), sep = ".")
    ggsave(
      file.path(obj@params$out_path, paste(prefixName, "sankey", ifelse(treshold > 0, "th", ""), "png", sep = ".")),
      pl,
      width = 8,
      height = 6
    )

    # HEATMAP
    # remove from table all datapoints with Archetype and also factors (i dont wont row and columns to 0)
    # df <- table(data$type[-"Archetype"], data$archetype[-"Archetype"])
    df <- table(data$type, data$archetype)
    df <- df[, colnames(df) != "Archetype"]
    df <- df[row.names(df) != "Archetype", ]
    df <- as.data.frame(df)

    # Heatmap with text inside
    plt_hm <- ggplot(df, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = Freq), color = "white") +
      geom_text(aes(label = Freq), vjust = 1) +
      scale_fill_gradient(low = "white", high = "blue") +
      theme_alluvial() +
      labs(x = "Cell types", y = "Archetype") +
      theme(axis.text.x = element_text(hjust = 1))
    plt_hm

    ggsave(
      file.path(obj@params$out_path, paste(prefixName, "heatmap", ifelse(treshold > 0, "th", ""), "png", sep = ".")),
      width = 8,
      height = 6
    )
  }
  # end treshold 2

  # Heatmap
}
# end k

