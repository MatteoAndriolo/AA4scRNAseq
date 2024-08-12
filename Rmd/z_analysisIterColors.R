library(Rtsne)
library(plotly)
library(ggplot2)
library(cowplot)
if (!require(ggsankey)) devtools::install_github("davidsjoberg/ggsankey")
if (!require(networkD3)) install.packages("networkD3")
library(ggsankey)
library(viridis)
library(Seurat)
library(scales)
library(dplyr)
if (!require(tidyverse)) install.packages("tidyverse")

set.seed(2024)

se <- readRDS(file.path("data/Melanoma/se.Rds"))
name <- list(
  # "GLYK",
  # "MAPK",
  # "CANCER",
  # "MTOR",
  # "TGF"
  "HVF"
)

pathw <- 1
treshold <- FALSE
out_path_base <- "out/Melanoma/images"
in_path <- "data/Melanoma"

# Function Definitions ---------------------------------------------------------
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
  colnames(weights) <- paste("Archetype", 1:num_archetypes, sep = "_")

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
  red <- "umap"
  treshold <- FALSE
}
for (k in c("7", "12")) {
  message("Using k=", k)
  num_archetypes <- as.integer(k)
  for (i in 1:length(name)) {
    pathw <- i
    namePathw <- name[pathw]
    message("Pathw ", namePathw, " with k ", k)
    out_path <- file.path(out_path_base, namePathw)
    if (!dir.exists(out_path)) dir.create(out_path)

    se3D <- readRDS(file.path(in_path, paste0(namePathw, ".Rds")))
    aa <- readRDS(file.path(in_path, paste0(namePathw, "_AA.Rds")))

    ################################################################################
    # newse <- addArchetypesToSeurat(se, aa, "7")
    newse <- findClosestPoints(se, se3D, aa, k)
    ################################################################################

    # colors_types <- rainbow(length(unique(newse$ctype)) - 1)
    num_types <- length(unique(newse$ctype))
    all_colors <- rainbow(num_types + num_archetypes - 1)
    colors_types <- c(all_colors[1:(num_types - 1)], "black")
    colors_archetypes <- all_colors[num_types:(num_types + num_archetypes - 1)]
    colors_Archetype <- "black"
    colors_text_Archetype <- c("white")

    size_types <- 1
    size_Archetype <- 4
    size_text_Archetype <- 3


    for (red in c("umap", "tsne")) {
      plot_data <- as.data.frame(Embeddings(newse, reduction = red))
      colnames(plot_data) <- c("X1", "X2", "X3")
      plot_data$Label <- newse$ctype

      message("reduction ", red)
      # 3D Plot
      p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = Label)) +
        geom_point(data = subset(plot_data, Label != "Archetype"), size = size_types) +
        geom_point(data = subset(plot_data, Label == "Archetype"), size = size_Archetype) +
        geom_text(
          data = subset(plot_data, Label == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = c(colors_types, colors_Archetype)) +
        theme_classic()
      p1

      p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = Label)) +
        geom_point(data = subset(plot_data, Label != "Archetype"), size = size_types) +
        geom_point(data = subset(plot_data, Label == "Archetype"), size = size_Archetype) +
        geom_text(
          data = subset(plot_data, Label == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = c(colors_types, colors_Archetype)) +
        theme_classic()
      p2

      p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = Label)) +
        geom_point(data = subset(plot_data, Label != "Archetype"), size = size_types) +
        geom_point(data = subset(plot_data, Label == "Archetype"), size = size_Archetype) +
        geom_text(
          data = subset(plot_data, Label == "Archetype"),
          aes(label = seq(1, num_archetypes)),
          color = colors_text_Archetype, size = size_text_Archetype
        ) +
        scale_color_manual(values = c(colors_types, colors_Archetype)) +
        theme_classic()
      p3

      combined_plot <- plot_grid(p1, p2, p3, ncol = 1)

      prefixName <- paste(name[pathw], k, red, sep = ".")
      ggsave(
        file.path(out_path, paste(prefixName, ".ct.X1.vs.X2.png")),
        p1,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(out_path, paste(prefixName, ".ct.X1.vs.X3.png")),
        p2,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(out_path, paste(prefixName, ".ct.X2.vs.X3.png")),
        plot = p3,
        width = 8,
        height = 6
      )
      ggsave(
        file.path(out_path, paste(prefixName, ".ct.png")),
        combined_plot,
        width = 8,
        height = 18
      )


      for (treshold in c(TRUE, FALSE)) {
        plot_data$aa_clusters <- as.character(newse@misc$aa_clusters)
        mycolors <- c(colors_types[1:length(unique(plot_data$aa_clusters))], "black")
        if (treshold) {
          plot_data$aa_clusters <- as.character(newse@misc$aa_cluster_treshold.5)
          mycolors <- c(colors_types[1:length(unique(plot_data$aa_clusters)) - 1], "black", "grey")
        }
        plot_data$aa_clusters[newse$ctype == "Archetype"] <- "Archetype"
        # aa clusters
        p1 <- ggplot(plot_data, aes(x = X1, y = X2, color = aa_clusters)) +
          geom_point(data = subset(plot_data, aa_clusters != "Archetype"), size = size_types) +
          geom_point(data = subset(plot_data, aa_clusters == "Archetype"), size = size_Archetype) +
          scale_color_manual(values = mycolors) +
          geom_text(
            data = subset(plot_data, aa_clusters == "Archetype"),
            aes(label = seq(1, num_archetypes)),
            color = colors_text_Archetype, size = size_text_Archetype
          ) +
          theme_classic()
        p1

        p2 <- ggplot(plot_data, aes(x = X1, y = X3, color = aa_clusters)) +
          geom_point(data = subset(plot_data, aa_clusters != "Archetype"), size = size_types) +
          geom_point(data = subset(plot_data, aa_clusters == "Archetype"), size = size_Archetype) +
          scale_color_manual(values = mycolors) +
          geom_text(
            data = subset(plot_data, aa_clusters == "Archetype"),
            aes(label = seq(1, num_archetypes)),
            color = colors_text_Archetype, size = size_text_Archetype
          ) +
          theme_classic()
        p2

        p3 <- ggplot(plot_data, aes(x = X2, y = X3, color = aa_clusters)) +
          geom_point(data = subset(plot_data, aa_clusters != "Archetype"), size = size_types) +
          geom_point(data = subset(plot_data, aa_clusters == "Archetype"), size = size_Archetype) +
          scale_color_manual(values = mycolors) +
          geom_text(
            data = subset(plot_data, aa_clusters == "Archetype"),
            aes(label = seq(1, num_archetypes)),
            color = colors_text_Archetype, size = size_text_Archetype
          ) +
          theme_classic()
        p3

        combined_plot <- plot_grid(p1, p2, p3, ncol = 1)

        prefixName <- paste(name[pathw], k, red, ifelse(treshold > 0, "th", ""), sep = ".")
        ggsave(
          file.path(out_path, paste(prefixName, ".aa.X1.vs.X2.png")),
          p1,
          width = 8,
          height = 6
        )
        ggsave(
          file.path(out_path, paste(prefixName, ".aa.X1.vs.X3.png")),
          p2,
          width = 8,
          height = 6
        )
        ggsave(
          file.path(out_path, paste(prefixName, ".aa.X2.vs.X3.png")),
          plot = p3,
          width = 8,
          height = 6
        )
        ggsave(
          file.path(out_path, paste(prefixName, ".aa.comb.png")),
          combined_plot,
          width = 8,
          height = 18
        )
      } # end treshold
    } # end red

    # 3D Plot
    p3d <- plot_ly(
      plot_data,
      x = ~X1, y = ~X2, z = ~X3,
      color = ~Label,
      colors = colors_types,
      size = ~ ifelse(Label == "Archetype", size_Archetype, size_types),
      type = "scatter3d",
      text = ~ ifelse(Label == "Archetype", seq(1, num_archetypes), ""),
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
    p3d

    # SANKEY PLOT
    # Assuming newse$ctype and newse@misc$aa_clusters are vectors of the same length
    for (treshold in c(TRUE, FALSE)) {
      ctype <- plot_data$Label
      taa_clusters <- as.character(newse@misc$aa_clusters)
      if (treshold) taa_clusters <- as.character(newse@misc$aa_cluster_treshold.5)

      # Create a data frame with the required columns
      data <- as.data.frame(list(
        type = plot_data$Label,
        archetype = taa_clusters
      ))
      d <- data %>%
        make_long(colnames(data)) %>%
        filter(node != "Archetype")

      pl <- ggplot(d, aes(
        x = x,
        next_x = next_x,
        node = node,
        next_node = next_node,
        fill = factor(node),
        label = node
      )) +
        geom_sankey(
          flow.alpha = 0.5, # This Creates the transparency of your node
          node.color = "black", # This is your node color
          show.legend = FALSE
        ) + # This determines if you want your legend to show
        geom_sankey_label(size = 3, color = "black", fill = "white") +
        labs(title = paste0("Sankey Tirosh - ", namePathw), subtitle = paste(ifelse(treshold, "With", "Without"), "treshold"), sep = " ") +
        theme_classic()


      pl
      prefixName <- paste(name[pathw], k, ifelse(treshold > 0, "th", ""), sep = ".")
      ggsave(
        file.path(out_path, paste0(prefixName, ".sankey.", ifelse(treshold > 0, ".th", ""), ".png")),
        pl,
        width = 8,
        height = 6
      )

      # HEATMAP
      df <- table(ctype, taa_clusters)
      df

      df_melt <- reshape2::melt(df, id.vars = "aa_clusters", variable.name = "ctype", value.name = "count")

      plt_hm <- ggplot(df_melt, aes(x = taa_clusters, y = ctype, fill = count)) +
        geom_tile(color = "white") +
        geom_text(aes(label = count), color = "black", size = 3) +
        scale_fill_gradient(low = "white", high = "blue") +
        theme_classic() +
        labs(
          title = "Heatmap of Cell Type Distribution Across Clusters", subtitle = paste(namePathw, "Tirosh", sep = " "),
          x = "Archetype", y = "Cell Type"
        ) +
        theme(axis.text.x = element_text(hjust = 1))
      plt_hm

      ggsave(
        file.path(out_path, paste0(prefixName, ".heatmap.", ifelse(treshold > 0, ".th", ""), ".png")),
        plt_hm,
        width = 8,
        height = 6
      )
    } # end treshold 2

    # Heatmap
  } # end i pathways
} # end k
