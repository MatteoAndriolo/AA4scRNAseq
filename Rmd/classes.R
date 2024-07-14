# Define the 'database' Class ----
setClass("database",
  slots = list(
    se = "Seurat",
    se.org = "ANY",
    plots = "list", # Group plots here
    archetypes = "list", # Group archetype analysis related stuff here
    params = "list", # Execution parameters
    curr.params = "list", # Current execution parameters
    compare = "list",
    other = "list"
  )
)
# DATA ----
## generics ------

### loadData ----
# Generic method for loading data
setGeneric("obj_loadData", function(obj,
                                    data_path = NULL,
                                    ...) {
  standardGeneric("obj_loadData")
})

### createSeuratObject ----
setGeneric("obj_createSeuratObject", function(obj, se, gene_names, cell_metadata, where.cell_names, pathw = NULL, test = FALSE, HVF = FALSE, test_genes = 300, test_samples = 500) {
  standardGeneric("obj_createSeuratObject")
})

### getSeData ----
setGeneric("obj_getSeData", function(obj) {
  standardGeneric("obj_getSeData")
})

### getMatrixHVF ----
setGeneric("obj_getMatrixHVF", function(obj) {
  standardGeneric("obj_getMatrixHVF")
})

# ### getCellTypesList ----
# setGeneric("obj_getCellTypesList", function(obj) {
#   standardGeneric("obj_getCellTypesList")
# })
#
# ### getCellTypesMetaDataName ----
# setGeneric("obj_getCellTypesMetaDataName", function(obj) {
#   standardGeneric("obj_getCellTypesMetaDataName")
# })

## updateParams ----
setGeneric("obj_updateParams", function(obj, updateCurrent = FALSE, ...) {
  standardGeneric("obj_updateParams")
})

setMethod("obj_updateParams", "database", function(obj, updateCurrent = FALSE, ...) {
  if (debug) message("DEBUG: obj_updateParams | entered")
  # for each key=value in ... do obj@params$key <- value
  list.params <- list(...)
  list.params
  for (i in names(list.params)) {
    message("DEBUG: obj_updateParams | Updating ", i, " with ", list.params[[i]])
    if (i == "pathw" & is.numeric(list.params[[i]])) {
      if (list.params[[i]] > 0) {
        if (debug) message("DEGBUG: obj_updateParams | Entering in is numeric with pathw ", list.params[[i]], " >0 ")
        obj@params$pathw <- pathways[[list.params[[i]]]]
      } else {
        if (debug) message("DEGBUG: obj_updateParams | Entering in is numeric with pathw ", list.params[[i]], "<=0")
        obj@params$pathw <- NULL
      }
    }
    obj@params[[i]] <- list.params[[i]]
  }

  return(obj)
})

## setGenes ----
setGeneric("obj_setGenes", function(obj, pathGenes = "/app/data/list_genes_pathway.RData") {
  standardGeneric("obj_setGenes")
})

setMethod("obj_setGenes", "database", function(obj, pathGenes = "/app/data/list_genes_pathway.RData") {
  if (debug) message("DEBUG: obj_setGenes | path genes is :", pathGenes)
  load(pathGenes)
  list_genes_human_pathway <- list_genes_human_pathway
  list_genes_mouse_pathway <- list_genes_mouse_pathway

  message("LOG: obj_setGenes | Setting Genes ", obj@params$pathw)

  if (inherits(obj, "Melanoma")) {
    #if (is.character(obj@params$pathw) && length(obj@params$pathw) == 1) {
      # If only one pathwas
      genes <- list_genes_human_pathway[[obj@params$pathw]]
      if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", length(genes))
    #} else if (is.list(obj@params$pathw) || is.vector(obj@params$pathw)) {
    #  # if more than one
    #  genes <- lapply(obj@params$pathw, function(p) list_genes_human_pathway[[p]])
    #  if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", sapply(genes, length))
    #}
  } else {
    #if (is.character(obj@params$pathw) && length(obj@params$pathw) == 1) {
      # If only one pathwas
      genes <- list_genes_mouse_pathway[[obj@params$pathw]][[1]]
      if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", length(genes))
    #} else if (is.list(obj@params$pathw) || is.vector(obj@params$pathw)) {
    #  # if more than one
    #  genes <- lapply(obj@params$pathw, function(p) list_genes_mouse_pathway[[p]])
    #  if (debug) message("DEBUG: obj_setGenes | Number genes in ", pathw, " is ", sapply(genes, length))
    #}
  }

  if (is.null(genes)) {
    error("Pathway not found")
    stop("Pathway not found")
  }

  obj <- obj_updateParams(obj,
    updateCurrent = TRUE,
    genes = genes
  )
  return(obj)
})

# Visualization ----

## visualizeData ----
# Method to visualize data
setGeneric("obj_visualizeData", function(obj) {
  standardGeneric("obj_visualizeData")
})

setMethod("obj_visualizeData", "database", function(obj) {
  path_figures <- obj@params$path_figures

  # PCA UMAP
  obj@plots$pca <- PCAPlot(obj@se)
  obj@plots$umap <- UMAPPlot(obj@se)

  obj@plots$combined_plot <- plot_grid(
    obj@plots$pca + theme(legend.position = "none"),
    obj@plots$umap + theme(legend.position = "none"),
    labels = c("A", "B"),
    ncol = 2
  )

  ggsave(filename = file.path(path_figures, "PCAPlot.png"), plot = obj@plots$pca)
  ggsave(filename = file.path(path_figures, "UMAPPlot.png"), plot = obj@plots$umap)
  ggplot2::ggsave(filename = file.path(path_figures, "combined_plot.png"), plot = obj@plots$combined_plot)

  # ELBOWPLOT
  obj@plots$elbowplot <- ElbowPlot(obj@se)
  ggsave(filename = file.path(path_figures, "ElbowPlot.png"), plot = obj@plots$elbowplot)

  # GOLD
  # setMethod("obj_plotGoldUmap", "Melanoma", function(obj) {
  #   ct <- "non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK."
  #   umap_celltypes <- DimPlot(obj@se, reduction = "umap", group.by = ct)
  #   umap_tumor <- DimPlot(obj@se, reduction = "umap", group.by = "tumor")
  #   return(list(umap_celltypes=umap_celltypes, umap_tumor=umap_tumor))
  # })
  newplots <- obj_plotGoldUmap(obj)
  do.call(function(...) {
    for (plot_name in names(list(...))) {
      obj@plots[[plot_name]] <- list(...)[[plot_name]]
      ggsave(filename = file.path(path_figures, paste0(plot_name, ".png")), plot = obj@plots[[plot_name]])
    }
  }, newplots)

  obj@plots$umap_seucl <- DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
  obj@plots$umap_aacl <- DimPlot(obj@se, reduction = "umap", group.by = "aa_clusters")

  return(obj)
})


## visualizeArchetypes ----
setGeneric("obj_visualizeArchetypes", function(obj) {
  standardGeneric("obj_visualizeArchetypes")
})

setMethod("obj_visualizeArchetypes", "database", function(obj) {
  # a <- obj@archetypes$model
  # k <- a$k
  path_figures <- obj@params$path_figures
  plotarchetyps <- xyplot(obj@archetypes$model, as.matrix(obj_getSeData(obj)))
  obj@plots$xyplot <- plotarchetyps
  ggsave(filename = file.path(path_figures, "archetypes_xyplot.png"), plot = obj@plots$xyplot)
  return(obj)
})

## umapArchetypes ----
setGeneric("obj_umapArchetypes", function(obj, treshold = 0.1) {
  standardGeneric("obj_umapArchetypes")
})

setMethod("obj_umapArchetypes", "database", function(obj, treshold = 0.01) {
  if (debug) message("DEBUG: obj_umapArchetypes | entering function ")
  umap_result <- UMAPPlot(obj@se)
  umap_data <- as.data.frame(umap_result$data)[, 1:2]
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  if (debug) message("DEBUG: obj_umapArchetypes | umapPlot done and fetched data")

  weights <- coef(obj@archetypes$model)
  weights <- as.data.frame(weights)
  if (debug) message("DEBUG: obj_umapArchetypes | weights dimension is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  weights[weights < treshold] <- 0

  # column_sums <- colSums(obj@archetypes$model$archetypes)
  # if (debug) message("DEBUG: obj_umapArchetypes | column_sums dimension is ", dim(column_sums)[[1]], " ", dim(column_sums)[[2]])
  # normalized_mat <- sweep(obj@archetypes$model$archetypes, 2, column_sums, FUN = "/")
  # if (debug) message("DEBUG: obj_umapArchetypes | normalized_mat dimension is ", dim(normalized_mat)[[1]], " ", dim(normalized_mat)[[2]])
  # weights <- as.data.frame(normalized_mat)

  if (debug) message("DEBUG: obj_umapArchetypes | weights dimension is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  if (debug) message("DEBUG: obj_umapArchetypes | umap_data dimension is ", dim(umap_data)[[1]], " ", dim(umap_data)[[2]])

  plot_list <- list()

  for (i in 1:obj@params$k) {
    if (debug) message("DEBUG: obj_umapArchetypes | weights dimension is ", length(weights[, i]))
    umap_data$weight <- weights[, i]
    plot_title <- sprintf("UMAP Archetype %d", i)
    umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = weight)) +
      geom_point(size = 1) +
      scale_color_gradient(low = "grey", high = "red") +
      ggtitle(plot_title) +
      labs(color = "Weight")
    # print(umap_plot)

    plot_list[[i]] <- umap_plot
  }

  combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
  obj@plots$umap_archetypes <- combined_plot

  obj <- obj_updateParams(obj, updateCurrent = TRUE, umap_threshold = treshold)
  return(obj)
})

## umapWithArchetypes ----
setGeneric("obj_umapWithArchetypes", function(obj, treshold = 0.1) {
  standardGeneric("obj_umapWithArchetypes")
})

setMethod("obj_umapWithArchetypes", "database", function(obj, treshold = 0.01) {
  message("LOG: obj_umapWithArchetypes | creating plot")
  if (debug) message("DEBUG: obj_umapWithArchetypes | entering function ")
  if (debug) message("DEBUG: obj_umapWithArchetypes | archetypes dimension is ", dim(parameters(obj@archetypes$model))[[1]], " ", dim(parameters(obj@archetypes$model))[[2]])
  aspe <- t(parameters(obj@archetypes$model))
  rownames(aspe) <- rownames(obj@se@assays$RNA$counts)
  colnames(aspe) <- paste0("Archetype", 1:ncol(aspe))

  # Combine the original matrix and archetypes
  newse <- cbind(as.matrix(obj@se@assays$RNA$counts), aspe)
  newse <- as(newse, "dgCMatrix")
  if (debug) message("DEBUG: obj_umapWithArchetypes | newse dimension is ", dim(newse)[[1]], " ", dim(newse)[[2]])

  # Create a temporary Seurat object with the combined matrix
  combined_obj <- CreateSeuratObject(counts = newse)
  # TODO check if this must be done!!!
  # combined_obj <- ScaleData(combined_obj, layer = "counts")
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

  obj@plots$umap_with_archetypes <- plot2
  message("LOG: obj_umapWithArchetypes | finished plot")
  return(obj)
})

### plotGoldUmap - generic ----
setGeneric("obj_plotGoldUmap", function(obj) {
  standardGeneric("obj_plotGoldUmap")
})

# General ----
## nameFiles -----
setGeneric("obj_nameFiles", function(obj, name, ext) {
  standardGeneric("obj_nameFiles")
})

setMethod("obj_nameFiles", "database", function(obj, name, ext) {
  if (obj@params$test) {
    t <- "T"
  } else {
    t <- ""
  }

  if (obj@params$hvf) {
    h <- "H"
  } else {
    h <- ""
  }
  if (is.null(obj@params$pathw)) {
    p <- ""
  } else {
    p <- sprintf("%s%s", "_", substr(obj@params$pathw, 1, 4))
  }

  return(sprintf("%s/%s%s%s%s_%s.%s", obj@params$out_path, class(obj), p, t, h, name, ext))
})

## seuratCluster ----
setGeneric("obj_seuratCluster", function(obj) {
  standardGeneric("obj_seuratCluster")
})

setMethod("obj_seuratCluster", "database", function(obj) {
  obj@se <- Seurat::FindNeighbors(obj@se, dims = 1:10)
  obj@se <- Seurat::FindClusters(obj@se, method = "igraph", resolution = 1, n.start = 10, n.iter = 10, verbose = TRUE)
  return(obj)
})


## saveObj ----
# Method to save object
setGeneric("obj_saveObj", function(obj, namefile = "final", keep.org = FALSE) {
  standardGeneric("obj_saveObj")
})

setMethod("obj_saveObj", "database", function(obj, namefile = "", keep.org = FALSE) {
  # filename <- sprintf("%s/%s_%s.rds", obj@params$out_path, class(obj), substr(obj@params$pathw, 1, 4))
  filename <- obj_nameFiles(obj, namefile, "rds")
  message(sprintf("LOG: obj_saveObj | Saving object to %s", filename))
  t <- obj
  if (!keep.org) {
    t@se.org <- NULL
  }
  t@curr.params <- list()
  saveRDS(t, file = filename)
})


# Other sources ----
source("/app/Rmd/class_Melanoma.R")
source("/app/Rmd/class_Exp.R")
source("/app/Rmd/class_archetypes.R")
# source("/app/Rmd/class_archetypal.R")
# source("/app/Rmd/class_Other.R")
