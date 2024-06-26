## General ----
### obj_loadData ----
# Define the 'database' Class
setClass("database",
  slots = list(
    se = "Seurat",
    se.org = "ANY",
    # data = "list", # Group data matrices here
    plots = "list", # Group plots here
    archetypes = "list", # Group archetype analysis related stuff here
    params = "list", # Execution parameters
    curr.params = "list", # Current execution parameters
    compare = "list"
  )
)

### obj_areParamsEqual  ----
# all but pathw that is ok if it changes
# setGeneric("obj_areParamsEqual", function(obj, update = FALSE) {
#  standardGeneric("obj_areParamsEqual")
# })
#
# setMethod("obj_areParamsEqual", "database", function(obj, update = FALSE) {
#  # check if all parameters (but pathw) are the same
#  flag <- TRUE
#  for (i in names(obj@params)) {
#    if (i != "pathw") {
#      if (!(is.null(obj@params[[i]]) | is.null(obj@curr.params[[i]])) & obj@params[[i]] != obj@curr.params[[i]]) {
#        flag <- FALSE
#        if (update) {
#          obj@curr.params[[i]] <- obj@params[[i]]
#        }
#      }
#    }
#  }
#  return(flag)
# })

### obj_updateParams ----
setGeneric("obj_updateParams", function(obj, updateCurrent = FALSE, ...) {
  standardGeneric("obj_updateParams")
})

setMethod("obj_updateParams", "database", function(obj, updateCurrent = FALSE, ...) {
  message("LOG: obj_updateParams | entered")
  # for each key=value in ... do obj@params$key <- value
  list.params <- list(...)
  list.params
  for (i in names(list.params)) {
    message("DEBUG: obj_updateParams | Updating ", i, " with ", list.params[[i]])
    if (i == "pathw" & is.numeric(list.params[[i]])) {
      if (list.params[[i]] > 0) {
        if (debug) warning("DEGBUG: obj_updateParams | Entering in is numeric with pathw ", list.params[[i]], ">0")
        obj@params$pathw <- pathways[[obj@params$pathw]]
      } else {
        if (debug) warning("DEGBUG: obj_updateParams | Entering in is numeric with pathw ", list.params[[i]], "<=0")
        obj@params$pathw <- NULL
      }
    }
    obj@params[[i]] <- list.params[[i]]
  }

  # if (!is.logical(updateCurrent)) {
  #  message("updateCurrent must be a logical value but it is ", updateCurrent, " setting it to FALSE")
  #  updateCurrent <- FALSE
  # }
  # if (updateCurrent) {
  #  obj@curr.params <- obj@params
  # }
  return(obj)
})

### obj_loadData ----
# Generic method for loading data
setGeneric("obj_loadData", function(obj,
                                    data_path = NULL,
                                    ...) {
  standardGeneric("obj_loadData")
})

### obj_createSeuratObject
setGeneric("obj_createSeuratObject", function(obj, se, gene_names, cell_metadata, where.cell_names, pathw = NULL, test = FALSE, HVF = FALSE, test_genes = 300, test_samples = 500) {
  standardGeneric("obj_createSeuratObject")
})

### obj_getSeData ----
setGeneric("obj_getSeData", function(obj) {
  standardGeneric("obj_getSeData")
})

### obj_setGenes ----
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
    if (is.character(obj@params$pathw) && length(obj@params$pathw) == 1) {
      genes <- list_genes_human_pathway[[obj@params$pathw]]
      if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", length(genes))
    } else if (is.list(obj@params$pathw) || is.vector(obj@params$pathw)) {
      genes <- lapply(obj@params$pathw, function(p) list_genes_human_pathway[[p]])
      if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", sapply(genes, length))
    }
  } else {
    if (is.character(obj@params$pathw) && length(obj@params$pathw) == 1) {
      genes <- list_genes_mouse_pathway[[obj@params$pathw]]
      if (debug) message("DEBUG: obj_setGenes | Number genes in ", obj@params$pathw, " is ", length(genes))
    } else if (is.list(obj@params$pathw) || is.vector(obj@params$pathw)) {
      genes <- lapply(obj@params$pathw, function(p) list_genes_mouse_pathway[[p]])
      if (debug) message("DEBUG: obj_setGenes | Number genes in ", pathw, " is ", sapply(genes, length))
    }
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

### obj_getMatrixHVF ----
setGeneric("obj_getMatrixHVF", function(obj) {
  standardGeneric("obj_getMatrixHVF")
})

### obj_visualizeData ----
# Method to visualize data
setGeneric("obj_visualizeData", function(obj) {
  standardGeneric("obj_visualizeData")
})

setMethod("obj_visualizeData", "database", function(obj) {
  obj@plots$pca <- PCAPlot(obj@se)

  obj@plots$umap <- UMAPPlot(obj@se)

  obj@plots$combined_plot <- plot_grid(
    obj@plots$pcaplot + theme(legend.position = "none"),
    obj@plots$umapplot + theme(legend.position = "none"),
    labels = c("A", "B"),
    ncol = 2
  )
  obj@plots$elbowplot <- ElbowPlot(obj@se)

  # print(combined_plot)
  # print(elbowplot)

  return(obj)
})

### obj_furthestSum ----
setGeneric("obj_furthestSum", function(obj, k = 5) {
  standardGeneric("obj_furthestSum")
})

# no obj@data
setMethod("obj_furthestSum", "database", function(obj, k = 5) {
  irows <- archetypal::find_furthestsum_points(obj@data$m, k = k, nworkers = parallell::detectCores() / 2)
  return(irows)
})

### obj_performArchetypes ----
# Method to perform archetypal analysis
setGeneric("obj_performArchetypes", function(obj, k = NULL, doparallel = TRUE) {
  standardGeneric("obj_performArchetypes")
})

setMethod("obj_performArchetypes", "database", function(obj, k = NULL, doparallel = FALSE) {
  message("LOG: obj_performArchetyps | Performing Archetypes and ", obj@params$pathw)
  if(k=NULL){
    k <- kneedle(obj@plots$elbowplot$data$dims, obj@plots$elbowplot$data$stdev)[1]
    message("LOG: obj_performArchetypes | Number of archetypes is ", k)
  }

  m <- as.matrix(obj_getSeData(obj))
  m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
  message("LOG: obj_performArchetyps | MATRIX DIMENSION FOR ARCHETYPES IS ", dim(m)[[1]], " ", dim(m)[[2]])

  #  obj@data$m <- as.matrix(obj_getSeData(obj))
  #  obj@data$m <- obj@data$m[Matrix::rowSums(obj@data$m) > 0, Matrix::colSums(obj@data$m) > 0]
  #  message("LOG: MATRIX DIMENSION FOR ARCHETYPES IS ", dim(obj@data$m)[[1]], " ", dim(obj@data$m)[[2]])
  #### Furthest Sum Initialization
  #  message("LOG: Finding Furthest Sum")
  #  irows <- obj_furthestSum(obj, k = k)
  #  message("LOG: Furthest Sum Found")

  #### Archetypes Computation
  obj@archetypes$restarts <- list()
  obj <- obj_updateParams(obj,
    updateCurrent = TRUE,
    k = k
  )

  max_iterations <- obj@params$max_iterations
  num_restarts <- obj@params$num_restarts

  runArchetypes <- function(i, data, k, max_iterations) {
    message("LOG: obj_performArchetypes | Starting rerun ", i, "/", num_restarts)
    temp <- list()
    # tryCatch(
    #  {
    tstart <- Sys.time()
    family <- archetypes::archetypesFamily(which = "robust")
    temp$a <- archetypes::archetypes(data, k = k, verbose = TRUE, maxIterations = max_iterations, saveHistory = TRUE, family = family)
    tend <- Sys.time()
    message(sprintf("Archetypes Computed in %s", tend - tstart))

    temp$rss <- temp$a$rss
    temp$time <- tend - tstart
    #  },
    #  error = function(e) {
    #    temp$a <- NULL
    #    temp$rss <- Inf
    #    temp$time <- NA
    #    message(sprintf("Error in archetypes: %s", e$message))
    #  }
    # )
    return(temp)
  }

  tstartReruns <- Sys.time()
  if (doparallel) {
    nworkers <- parallel::detectCores() - 1
    # results <- mclapply(1:num_restarts, runArchetypes, data = obj@data$m, k = k, max_iterations = max_iterations, mc.cores = 3)
    results <- mclapply(1:num_restarts, runArchetypes, data = m, k = k, max_iterations = max_iterations, mc.cores = 3)
    obj@archetypes$restarts <- results
  } else {
    for (i in 1:num_restarts) {
      # obj@archetypes$restarts[[i]] <- runArchetypes(i, data = obj@data$m, k = k, max_iterations = max_iterations)
      message("LOG: obj_performArchetypes | Starting rerun ", i, "/", num_restarts)
      temp <- list()

      # tryCatch(
      # {
      tstart <- Sys.time()
      family <- archetypes::archetypesFamily(which = "robust") # , initfn = make.fix.initfn(irows[1]))
      # temp$a <- archetypes::archetypes(obj@data$m, k = k, verbose = TRUE, maxIterations = max_iterations, saveHistory = TRUE, family = family)
      temp$a <- archetypes::archetypes(m, k = k, verbose = TRUE, maxIterations = max_iterations, saveHistory = TRUE, family = family)
      tend <- Sys.time()
      message(sprintf("Archetypes Computed in %s", tend - tstart))
      temp$rss <- temp$a$rss
      temp$time <- tend - tstart
      # },
      # error = function(e) {
      #   temp$a <- NULL
      #   temp$rss <- Inf
      #   temp$time <- NA
      #   message(sprintf("Error in archetypes: %s", e$message))
      # }
      # )
      obj@archetypes$restarts[[i]] <- temp
    }
  }
  tendReruns <- Sys.time()

  message("LOG: obj_performArchetypes | Reruns completed in ", tendReruns - tstartReruns, " seconds")
  obj@archetypes$bestrun <- obj@archetypes$restarts[[which.min(sapply(obj@archetypes$restarts, function(x) x$rss))]]
  # obj@archetypes$restarts <- list()
  obj@archetypes$model <- obj@archetypes$bestrun$a

  return(obj)
})

### obj_visualizeArchetypes ----
setGeneric("obj_visualizeArchetypes", function(obj) {
  standardGeneric("obj_visualizeArchetypes")
})

setMethod("obj_visualizeArchetypes", "database", function(obj) {
  # a <- obj@archetypes$model
  # k <- a$k
  plotarchetyps <- xyplot(obj@archetypes$model, as.matrix(obj_getSeData(obj)))
  obj@plots$xyplot <- plotarchetyps
  return(obj)
})

### obj_umapArchetypes ----
setGeneric("obj_umapArchetypes", function(obj, treshold = 0.2) {
  standardGeneric("obj_umapArchetypes")
})

setMethod("obj_umapArchetypes", "database", function(obj, treshold = 0.2) {
  if (debug) message("DEBUG: obj_umapArchetypes | entering function ")
  umap_result <- UMAPPlot(obj@se)
  umap_data <- as.data.frame(umap_result$data)[, 1:2]
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  if (debug) message("DEBUG: obj_umapArchetypes | umapPlot done and fetched data")

  weights <- coef(obj@archetypes$model)
  weights <- as.data.frame(weights)
  if (debug) message("DEBUG: obj_umapArchetypes | weights dimension is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  weights[weights < treshold] <- 0


  column_sums <- colSums(obj@archetypes$model$archetypes)
  if (debug) message("DEBUG: obj_umapArchetypes | column_sums dimension is ", dim(column_sums)[[1]], " ", dim(column_sums)[[2]])
  normalized_mat <- sweep(obj@archetypes$model$archetypes, 2, column_sums, FUN = "/")
  if (debug) message("DEBUG: obj_umapArchetypes | normalized_mat dimension is ", dim(normalized_mat)[[1]], " ", dim(normalized_mat)[[2]])
  weights <- as.data.frame(normalized_mat)

  if (debug) message("DEBUG: obj_umapArchetypes | weights dimension is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  if (debug) message("DEBUG: obj_umapArchetypes | umap_data dimension is ", dim(umap_data)[[1]], " ", dim(umap_data)[[2]])

  plot_list <- list()
  
  for (i in 1:obj@params$k) {
    umap_data$weight <- t(weights[i, ])
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
  # print(combined_plot)
  obj@plots$umap_archetypes <- combined_plot

  obj <- obj_updateParams(obj, updateCurrent = TRUE, umap_threshold = treshold)
  # obj@params$umap_threshold <- treshold
  return(obj)
})

### obj_assignAAClusters ----
setGeneric("obj_assignAAClusters", function(obj) {
  standardGeneric("obj_assignAAClusters")
})

setMethod("obj_assignAAClusters", "database", function(obj) {
  # se <- obj@se
  # a <- obj@archetypes$model
  # k <- a$k
  message("LOG: obj_assingAACluster | creating aa_clusters metadata")
  # weights <- coef(obj@archetypes$model)
  weights <- obj@archetypes$model$archetypes
  if (debug) message("DEBUG: obj_assignAAClusters | dimension of weights is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  weights <- as.data.frame(weights)


  if (debug) message("LOG: obj_assignAAClusters | dimension of meta.data is ", dim(obj@se@meta.data)[[1]], " ", dim(obj@se@meta.data)[[2]])
  obj@se@meta.data$aa_clusters <- apply(weights, 2, which.max)
  message("LOG: obj_assingAACluster | finished aa_clusters metadata")
  return(obj)
})


### obj_nameFiles -----
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

  return(sprintf("%s/%s%s%s_%s.%s", obj@params$out_path, class(obj), t, h, name, ext))
})

### obj_seuratCluster ----
setGeneric("obj_seuratCluster", function(obj) {
  standardGeneric("obj_seuratCluster")
})

setMethod("obj_seuratCluster", "database", function(obj) {
  obj@se <- Seurat::FindNeighbors(obj@se, dims = 1:10)
  obj@se <- Seurat::FindClusters(obj@se, method = "igraph", resolution = 1, n.start = 10, n.iter = 10, verbose = TRUE)
  names(obj@se)
  obj@plots$clusterplot <- Seurat::DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
  return(obj)
})


### obj_saveObj ----
# Method to save object
setGeneric("obj_saveObj", function(obj, namefile = "final", keep.org = FALSE) {
  standardGeneric("obj_saveObj")
})

setMethod("obj_saveObj", "database", function(obj, namefile = "", keep.org = FALSE) {
  #filename <- sprintf("%s/%s_%s.rds", obj@params$out_path, class(obj), substr(obj@params$pathw, 1, 4))
  filename <- obj_nameFiles(obj, namefile, "rds")
  message(sprintf("LOG: obj_saveObj | Saving object to %s", filename))
  t <- obj
  if (!keep.org) {
    t@se.org <- NULL
  }
  t@curr.params <- list()
  saveRDS(t, file = filename)
})

source("/app/Rmd/class_Melanoma.R")
source("/app/Rmd/class_Exp.R")
# source("/app/Rmd/class_Other.R")
