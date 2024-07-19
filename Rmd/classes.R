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
    # if (is.character(obj@params$pathw) && length(obj@params$pathw) == 1) {

    # If only one pathwas
    genes <- list_genes_human_pathway[[obj@params$pathw]]
    if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", length(genes))

    # } else if (is.list(obj@params$pathw) || is.vector(obj@params$pathw)) {
    #  # if more than one
    #  genes <- lapply(obj@params$pathw, function(p) list_genes_human_pathway[[p]])
    #  if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", sapply(genes, length))
    # }
  } else {
    # if (is.character(obj@params$pathw) && length(obj@params$pathw) == 1) {

    # If only one pathwas
    genes <- list_genes_mouse_pathway[[obj@params$pathw]]
    if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", length(genes))

    # } else if (is.list(obj@params$pathw) || is.vector(obj@params$pathw)) {
    #  # if more than one
    #  genes <- lapply(obj@params$pathw, function(p) list_genes_mouse_pathway[[p]])
    #  if (debug) message("DEBUG: obj_setGenes | Number genes in ", pathw, " is ", sapply(genes, length))
    # }
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

#  obj@plots$umap_seucl <- DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
#  obj@plots$umap_aacl <- DimPlot(obj@se, reduction = "umap", group.by = "aa_clusters")

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
source("/app/Rmd/class_Mouse.R")
source("/app/Rmd/class_archetypes.R")
source("/app/Rmd/class_archetypal.R")
# source("/app/Rmd/class_Other.R")
