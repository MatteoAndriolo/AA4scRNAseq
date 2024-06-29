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
    compare = "list",
    other = "list"
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
    obj@plots$pca + theme(legend.position = "none"),
    obj@plots$umap + theme(legend.position = "none"),
    labels = c("A", "B"),
    ncol = 2
  )

  obj@plots$elbowplot <- ElbowPlot(obj@se)

  # obj@se@meta.data$seurat_clusters
  obj <- obj_plotObjSpecificUmap(obj)
  obj@plots$umap_tumor <- DimPlot(obj@se, reduction = "umap", group.by = "tumor")
  obj@plots$umap_seucl <- DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
  obj@plots$umap_aacl <- DimPlot(obj@se, reduction = "umap", group.by = "aa_clusters")

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
setGeneric("obj_performArchetypes", function(obj, kappas = NULL, k = NULL, doparallel = TRUE) {
  standardGeneric("obj_performArchetypes")
})

setMethod("obj_performArchetypes", "database", function(obj, kappas = NULL, k = NULL, doparallel = FALSE) {
  if (debug) {
    message("DEBUG: obj_performArchetypes | k=", k)
    message("DEBUG: obj_performArchetypes | kappas=", kappas)
    message("DEBUG: obj_performArchetypes | obj@params$k=", obj@params$k)
    message("DEBUG: obj_performArchetypes | obj@params$kappas=", obj@params$kappas)
  }
  if (is.null(obj@params$kappas) & is.null(obj@params$k)) {
    stop("ERROR: obj_performArchetypes | k and kappas are null")
  }
  if (is.null(obj@params$kappas) & !is.null(obj@params$k)) {
    obj@params$kappas <- obj@params$k
  }

  if (is.null(obj@params$which.aa)) {
    obj_updateParams(obj, which.aa = "robust")
  }

  # SETUP matrices
  message("LOG: obj_performArchetypes | Performing Archetypes on pathw ", obj@params$pathw)
  message("LOG: obj_performArchetypes | Number of archetypes is ", obj@params$kappas)

  # IMPORTANT !!!! column <-> features, row <-> samples
  m <- as.matrix(t(obj_getSeData(obj)))
  m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
  message("LOG: obj_performArchetypes | matrix dimension for archetypes is ", dim(m)[1], " ", dim(m)[2])


  #### Archetypes Computation
  obj@archetypes$aa.kappas <- list()

  runArchetypes <- function(i, data, k, max_iterations) {
    message("LOG: obj_performArchetypes | Starting rerun ", i, "/", num_restarts)
    temp <- list()
    # tryCatch(
    #  {
    tstart <- Sys.time()
    family <- archetypes::archetypesFamily(which = "robust")
    temp$a <- archetypes::archetypes(data, k = k, verbose = TRUE, maxIterations = max_iterations, saveHistory = FALSE, family = family)
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
  # if (doparallel) {
  #   nworkers <- parallel::detectCores() - 1
  #   # results <- mclapply(1:num_restarts, runArchetypes, data = obj@data$m, k = k, max_iterations = obj@params$max_iterations, mc.cores = 3)
  #   results <- mclapply(1:num_restarts, runArchetypes, data = m, k = k, max_iterations = obj@params$max_iterations, mc.cores = nworkers)
  #   obj@archetypes$restarts <- results
  # } else {
  family <- archetypes::archetypesFamily(which = obj@params$which.aa)
  obj@params$family <- family

  for (k in obj@params$kappas) {
    # Archetypal Analysis ---------------------------------------------------------
    history.restarts.k <- list()
    best_rss <- Inf
    # best_model=NULL
    best_restart_index <- -1

    for (i in 1:obj@params$num_restarts) {
      temp <- list()
      message("LOG: obj_performArchetypes | Starting rerun ", i, "/", obj@params$num_restarts, " with k=", k)
      # FURTHEST SUM ---------------------------------------------------------
      if (!is.null(obj@params$doFurthestSum) & obj@params$doFurthestSum) {
        message("LOG: obj_performArchetypes | Performing Furthest Sum")
        ttstart <- Sys.time()
        irow <- sample(1:nrow(m), 1)
        irows <- FurthestSum(m, irow = irow, k = k)
        family$initfn <- archetypes:::make.fix.initfn(irows)
        ttend <- Sys.time()
        temp$FStime <- difftime(ttend, ttstart, units = "secs")
        temp$FSindices <- sort(irows)

        message("LOG: obj_performArchetypes | Furthest Sum Done in ", temp$FStime, " seconds")
        
        message("LOG: obj_performArchetypes | Furthest Sum Done found ", paste(temp$FSindices, collapse = ", "))
      }


      # If you want to use function instead of explicit code uncomment this
      # obj@archetypes$restarts[[i]] <- runArchetypes(i, data = obj@data$m, k = k, max_iterations = obj@params$max_iterations)

      tstart <- Sys.time()
      temp$a <- archetypes::archetypes(m, k = k, verbose = TRUE, maxIterations = obj@params$max_iterations, saveHistory = TRUE, family = family)
      tend <- Sys.time()
      message("Archetypes Computed in ", difftime(tend, tstart, units = "secs"))

      temp$rss <- temp$a$rss
      temp$time <- difftime(tend, tstart, units = "secs")
      # },
      # error = function(e) {
      #   temp$a <- NULL
      #   temp$rss <- Inf
      #   temp$time <- NA
      #   message(sprintf("Error in archetypes: %s", e$message))
      # }
      # )
      history.restarts.k[[i]] <- temp

      # MISC
      message("DEBUG: obj_performArchetypes | best_rss is ", best_rss, " and temp$rss is ", temp$rss)
      if (is.na(temp$rss)) {
        temp$rss <- Inf
      }
      if (temp$rss < best_rss) {
        best_rss <- temp$rss
        # best_model = temp$a
        best_restart_index <- i
        if (debug) message("DEBUG: obj_performArchetypes | best_rss chosen is ", best_rss)
      }
    }
    history.restarts.k$best.run <- history.restarts.k[[best_restart_index]]
    obj@archetypes$aa.kappas[[as.character(k)]] <- history.restarts.k
  }
  # }
  tendReruns <- Sys.time()
  message("OUTPUT: obj_performArchetypes | Reruns completed in ", difftime(tendReruns, tstartReruns, units = "secs"), " seconds")

  # Now find the best model across all kappas
  best_overall_run <- NULL
  best_overall_rss <- Inf

  for (k in obj@params$kappas) {
    best_k_run <- obj@archetypes$aa.kappas[[as.character(k)]]$best.run
    if (best_k_run$rss < best_overall_rss) {
      best_overall_rss <- best_k_run$rss
      best_overall_run <- best_k_run
    }
  }
  # obj@archetypes$bestrun <- obj@archetypes$restarts[[which.min(sapply(obj@archetypes$restarts, function(x) x$rss))]]
  obj@archetypes$bestrun <- best_overall_run
  obj@archetypes$model <- best_overall_run$a
  if (debug) message("DEBUG: obj_performArchetypes | dim archetypes ", dim(parameters(best_overall_run$a))[1], " ", dim(parameters(best_overall_run$a))[2])
  # obj@archetypes$screeplot <- screeplot()
  # obj@archetypes$model <- obj@archetypes$bestrun$a

  # obj@archetypes$restarts <- list()

  return(obj)
})

# obj_performStepArchetypes
setGeneric("obj_performStepArchetypes", function(obj, kappas = NULL, k = NULL, doparallel = TRUE) {
  standardGeneric("obj_performStepArchetypes")
})

setMethod("obj_performStepArchetypes", "database", function(obj, kappas = NULL, k = NULL, doparallel = FALSE) {

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

### obj_umapWithArchetypes ----
setGeneric("obj_umapWithArchetypes", function(obj, treshold = 0.1) {
  standardGeneric("obj_umapWithArchetypes")
})

setMethod("obj_umapWithArchetypes", "database", function(obj, treshold = 0.01) {
  message("LOG: obj_umapWithArchetypes | creating plot")
  if (debug) message("DEBUG: obj_umapWithArchetypes | entering function ")
  if (debug) message("DEBUG: obj_umapWithArchetypes | archetypes dimension is ", dim(parameters(obj@archetypes$model))[[1]], " ", dim(parameters(obj@archetypes$model))[[2]])
  #################
  aspe <- t(parameters(obj@archetypes$model))
  rownames(aspe) <- rownames(obj@se@assays$RNA$counts)
  colnames(aspe) <- paste0("Archetype", 1:ncol(aspe))

  # Combine the original matrix and archetypes
  newse <- cbind(as.matrix(obj@se@assays$RNA$counts), aspe)
  newse <- as(newse, "dgCMatrix")
  if (debug) message("DEBUG: obj_umapWithArchetypes | newse dimension is ", dim(newse)[[1]], " ", dim(newse)[[2]])

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

  obj@plots$umap_with_archetypes <- plot2
  message("LOG: obj_umapWithArchetypes | finished plot")
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
  weights <- coef(obj@archetypes$model)
  if (debug) message("DEBUG: obj_assignAAClusters | dimension of weights is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  weights <- as.data.frame(weights)


  if (debug) message("LOG: obj_assignAAClusters | dimension of meta.data is ", dim(obj@se@meta.data)[[1]], " ", dim(obj@se@meta.data)[[2]])
  obj@se@meta.data$aa_clusters <- apply(weights, 1, which.max)
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
  if (is.null(obj@params$pathw)) {
    p <- ""
  } else {
    p <- sprintf("%s%s", "_", substr(obj@params$pathw, 1, 4))
  }

  return(sprintf("%s/%s%s%s%s_%s.%s", obj@params$out_path, class(obj), p, t, h, name, ext))
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

### obj_plotGoldUmap
setGeneric("obj_plotObjSpecificUmap", function(obj) {
  standardGeneric("obj_plotObjSpecificUmap")
})

### obj_getCellTypesList ----
setGeneric("obj_getCellTypesList", function(obj) {
  standardGeneric("obj_getCellTypesList")
})

### obj_getCellTypesMetaDataName ----
setGeneric("obj_getCellTypesMetaDataName", function(obj) {
  standardGeneric("obj_getCellTypesMetaDataName")
})

source("/app/Rmd/class_Melanoma.R")
source("/app/Rmd/class_Exp.R")
# source("/app/Rmd/class_Other.R")
