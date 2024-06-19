## General ----
### obj_loadData ----
# Define the 'database' Class
setClass("database",
  slots = list(
    se = "Seurat",
    se.org = "ANY",
    data = "list", # Group data matrices here
    plots = "list", # Group plots here
    archetypes = "list", # Group archetype analysis related stuff here
    params = "list", # Execution parameters
    curr.params = "list" # Current execution parameters
  )
)

### obj_setParams ----

# setGeneric("obj_setParams", function(obj, test = NULL, HVF = NULL, test_genes = NULL, test_samples = NULL, pathw = NULL, out_path = NULL, num_restarts = NULL, max_iterations = NULL) {
#  standardGeneric("obj_setParams")
# })
#
# setMethod("obj_setParams", "database", function(obj, test = NULL, HVF = NULL, test_genes = NULL, test_samples = NULL, pathw = NULL, out_path = NULL, num_restarts = NULL, max_iterations = NULL) {
#  if (!is.null(test)) {
#    obj@params$TEST<- test
#  }
#  if (!is.null(HVF)) {
#    obj@params$HVF <- HVF
#  }
#  if (!is.null(test_genes)) {
#    obj@params$test_genes <- test_genes
#  }
#  if (!is.null(test_samples)) {
#    obj@params$test_samples <- test_samples
#  }
#  if (!is.null(pathw)) {
#    obj@params$pathw <- pathw
#  }
#  if (!is.null(out_path)) {
#    obj@params$out_path <- out_path
#  }
#  if (!is.null(num_restarts)) {
#    obj@params$num_restarts <- num_restarts
#  }
#  if (!is.null(max_iterations)) {
#    obj@params$max_iterations <- max_iterations
#  }
#
#  if (is.null(obj@curr.params)) {
#    obj@curr.params <- obj@params
#  }
#  return(obj)
# })

### obj_areParamsEqual  ----
# all but pathw that is ok if it changes
setGeneric("obj_areParamsEqual", function(obj, update = FALSE) {
  standardGeneric("obj_areParamsEqual")
})

setMethod("obj_areParamsEqual", "database", function(obj, update = FALSE) {
  # check if all parameters (but pathw) are the same
  flag <- TRUE
  for (i in names(obj@params)) {
    if (i != "pathw") {
      if (!(is.null(obj@params[[i]]) | is.null(obj@curr.params[[i]])) & obj@params[[i]] != obj@curr.params[[i]]) {
        flag <- FALSE
        if (update) {
          obj@curr.params[[i]] <- obj@params[[i]]
        }
      }
    }
  }
  return(flag)
})

### obj_updateParams ----
setGeneric("obj_updateParams", function(obj, updateCurrent = FALSE, ...) {
  standardGeneric("obj_updateParams")
})

setMethod("obj_updateParams", "database", function(obj, updateCurrent = FALSE, ...) {
  # for each key=value in ... do obj@params$key <- value
  list.params <- list(...)
  for (i in names(list.params)) {
    if (debug) message("DEBUG: Updating ", i, " with ", list.params[[i]])
    obj@params[[i]] <- list.params[[i]]
  }

  if (!is.logical(updateCurrent)) {
    message("updateCurrent must be a logical value but it is ", updateCurrent, " setting it to FALSE")
    updateCurrent <- FALSE
  }
  if (updateCurrent) {
    obj@curr.params <- obj@params
  }
  return(obj)
})

### obj_loadData ----
# Generic method for loading data
setGeneric("obj_loadData", function(obj,
                                    data_path = NULL,
                                    test = FALSE,
                                    HVF = TRUE,
                                    pathw = NULL,
                                    test_genes = 300,
                                    test_samples = 500,
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
setGeneric("obj_setGenes", function(obj, pathw, pathGenes = "/app/data/list_genes_pathway.RData") {
  standardGeneric("obj_setGenes")
})

setMethod("obj_setGenes", "database", function(obj, pathw, pathGenes = "/app/data/list_genes_pathway.RData") {
  load(pathGenes)
  list_genes_human_pathway <- list_genes_human_pathway
  list_genes_mouse_pathway <- list_genes_mouse_pathway

  message("LOG: Setting Genes ", pathw)

  if (inherits(obj, "Melanoma")) {
    if (is.character(pathw) && length(pathw) == 1) {
      genes <- list_genes_human_pathway[[pathw]]
      if (debug) message("DEBUG: Number genes in ", pathw, " is ", length(genes))
    } else if (is.list(pathw) || is.vector(pathw)) {
      genes <- lapply(pathw, function(p) list_genes_human_pathway[[p]])
      if (debug) message("DEBUG: Number genes in ", pathw, " is ", sapply(genes, length))
    }
  } else {
    if (is.character(pathw) && length(pathw) == 1) {
      genes <- list_genes_mouse_pathway[[pathw]]
      if (debug) message("DEBUG: Number genes in ", pathw, " is ", length(genes))
    } else if (is.list(pathw) || is.vector(pathw)) {
      genes <- lapply(pathw, function(p) list_genes_mouse_pathway[[p]])
      if (debug) message("DEBUG: Number genes in ", pathw, " is ", sapply(genes, length))
    }
  }

  if (is.null(genes)) {
    error("Pathway not found")
    stop("Pathway not found")
  }

  obj <- obj_updateParams(obj,
    updateCurrent = TRUE,
    pathw = pathw,
    genes = genes,
    pathGenes = pathGenes
  )
  return(obj)
})

### obj_filterGenesFromOrg ----
setGeneric("obj_filterGenesFromOrg", function(obj, pathw) {
  standardGeneric("obj_filterGenesFromOrg")
})

setMethod("obj_filterGenesFromOrg", "database", function(obj, pathw) {
  message("Pathway")
  obj@se <- obj@se.org
  obj <- obj_setGenes(obj, pathw)
  message("LOG: Number of genes: ", length(obj@params$genes))
  obj@se <- obj@se[obj@params$genes, ]
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

setMethod("obj_furthestSum", "database", function(obj, k = 5) {
  irows <- archetypal::find_furthestsum_points(obj@data$m, k = k)
  return(irows)
})

### obj_performArchetypes ----
# Method to perform archetypal analysis
setGeneric("obj_performArchetypes", function(obj, k = NULL, HVF = FALSE, max_iters = 100, num_restarts = 1, doparallel=TRUE) {
  standardGeneric("obj_performArchetypes")
})

setMethod("obj_performArchetypes", "database", function(obj, k = NULL, HVF = FALSE, max_iters = 100, num_restarts = 1, doparallel=TRUE) {
  message("Performing Archetypes")
  if (HVF) {
    obj <- obj_getMatrixHVF(obj)
  }

  k <- kneedle(obj@plots$elbowplot$data$dims, obj@plots$elbowplot$data$stdev)[1]
  message("Number of archetypes is ", k)

  obj@data$m <- as.matrix(obj_getSeData(obj))
  obj@data$m <- obj@data$m[Matrix::rowSums(obj@data$m) > 0, Matrix::colSums(obj@data$m) > 0]

  #### Furthest Sum Initialization
  message("LOG: Finding Furthest Sum")
  irows <- obj_furthestSum(obj, k = k)
  message("LOG: Furthest Sum Found")

  #### Archetypes Computation

  obj@archetypes$restarts <- list()
  obj_updateParams(obj,
    updateCurrent = TRUE,
    k = k,
    max_iterations = max_iters,
    num_restarts = num_restarts
  )

  runArchetypes <- function(i, data, k, max_iterations) {
    message("Starting rerun ", i, "/", num_restarts)
    temp <- list()
    tryCatch(
      {
        tstart <- Sys.time()
        family <- archetypes::archetypesFamily(which = "robust")
        temp$a <- archetypes::archetypes(data, k = k, verbose = TRUE, maxIterations = max_iterations, saveHistory = TRUE, family = family)
        tend <- Sys.time()
        message(sprintf("Archetypes Computed in %s", tend - tstart))

        temp$rss <- temp$a$rss
        temp$time <- tend - tstart
      },
      error = function(e) {
        temp$a <- NULL
        temp$rss <- Inf
        temp$time <- NA
        message(sprintf("Error in archetypes: %s", e$message))
      }
    )
    return(temp)
  }



  tstartReruns <- Sys.time()
  if (doparallel) {
    nworkers=parallel::detectCores()-1
    results <- mclapply(1:num_restarts, runArchetypes, data = obj@data$m, k = k, max_iterations = max_iterations, mc.cores = nworkers)

    # Store the results in the object
    obj@archetypes$restarts <- results
  } else {
    for (i in 1:num_restarts) {
      runArchetypes(i, data = obj@data$m, k = k, max_iterations = max_iterations)
      #message("Starting rerun ", i, "/", num_restarts)
      #temp <- list()

      #tryCatch(
      #  {
      #    tstart <- Sys.time()
      #    family <- archetypes::archetypesFamily(which = "robust") # , initfn = make.fix.initfn(irows[1]))
      #    temp$a <- archetypes::archetypes(obj@data$m, k = k, verbose = TRUE, maxIterations = max_iterations, saveHistory = TRUE, family = family)
      #    tend <- Sys.time()
      #    message(sprintf("Archetypes Computed in %s", tend - tstart))
      #    temp$rss <- temp$a$rss
      #    temp$time <- tend - tstart
      #  },
      #  error = function(e) {
      #    temp$a <- NULL
      #    temp$rss <- Inf
      #    temp$time <- NA
      #    message(sprintf("Error in archetypes: %s", e$message))
      #  }
      #)
      #obj@archetypes$restarts[[i]] <- temp
    }
  }
  tendReruns <- Sys.time()

  message("Reruns completed in ", tendReruns - tstartReruns, " seconds")
  obj@archetypes$bestrun <- obj@archetypes$restarts[[which.min(sapply(obj@archetypes$restarts, function(x) x$rss))]]
  obj@archetypes$model <- obj@archetypes$bestrun$a

  return(obj)
})

### obj_visualizeArchetypes ----
setGeneric("obj_visualizeArchetypes", function(obj, out_path) {
  standardGeneric("obj_visualizeArchetypes")
})

setMethod("obj_visualizeArchetypes", "database", function(obj, out_path) {
  # a <- obj@archetypes$model
  # k <- a$k
  plotarchetyps <- xyplot(obj@archetypes$model, as.matrix(obj_getSeData(obj)))
  obj@plots$xyplot <- plotarchetyps
  return(obj)
})

### obj_umapArchetypes ----
setGeneric("obj_umapArchetypes", function(obj, out_path = NULL, treshold = 0.2) {
  standardGeneric("obj_umapArchetypes")
})

setMethod("obj_umapArchetypes", "database", function(obj, out_path = NULL, treshold = 0.2) {
  se <- obj@se
  a <- obj@archetypes$model
  k <- a$k

  umap_result <- UMAPPlot(se)
  umap_data <- as.data.frame(umap_result$data)[, 1:2]
  colnames(umap_data) <- c("UMAP1", "UMAP2")

  weights <- coef(a)
  weights <- as.data.frame(weights)
  weights[weights < treshold] <- 0

  column_sums <- colSums(a$archetypes)
  normalized_mat <- sweep(a$archetypes, 2, column_sums, FUN = "/")
  weights <- as.data.frame(normalized_mat)

  plot_list <- list()
  for (i in 1:k) {
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


### obj_saveObj ----
# Method to save object
setGeneric("obj_saveObj", function(obj) {
  standardGeneric("obj_saveObj")
})

setMethod("obj_saveObj", "database", function(obj) {
  filename <- sprintf("%s/%s/%s_%s.rds", obj@params$out_path, class(obj), class(obj), substr(obj@params$pathw, 1, 4))
  message(sprintf("Saving object to %s", filename))
  saveRDS(obj, file = filename)
})

# . #############################################################################
# Melanoma ----
# . #############################################################################
setClass("Melanoma",
  contains = "database"
)

### obj_loadData ----
setMethod(
  "obj_loadData", "Melanoma",
  function(obj,
           data_path = "/app/data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500,
           ...) {
    # obj <- obj_setParams(obj, test = test, HVF = HVF, test_genes = test_genes, test_samples = test_samples, pathw = pathw)
    # isnew <- obj_areParamsEqual(obj, update = TRUE)

    if (is.null(obj@se.org)) { #  isnew) {
      message("LOG: Full loading")
      se <- read.table(data_path, header = TRUE)
      se <- se[!duplicated(se[, 1]), ]
      rownames(se) <- se[, 1]
      se <- se[, 2:ncol(se)]

      metadata <- se[1:3, ]
      metadata <- t(metadata) %>%
        data.frame() %>%
        mutate(across(where(is.character), as.numeric))

      se <- se[4:nrow(se), ]
      se <- se %>%
        data.frame() %>%
        mutate(across(where(is.character), as.numeric))
      if (debug) message("DEBUG: data matrix has dimension ", dim(se)[[1]], " ", dim(se)[[2]])
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
      if (debug) message("DEBUG: data matrix has dimension post rem0", dim(se)[[1]], " ", dim(se)[[2]])


      if (test) {
        # tgenes <- min(test_genes, nrow(se))
        tgenes <- nrow(se)
        tsamples <- min(test_samples, ncol(se))
        metadata <- metadata[1:tsamples, ]

        se <- se[1:tgenes, 1:tsamples]
        se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
      }

      se <- CreateSeuratObject(counts = se, meta.data = metadata)
      message("LOG: Seurat object has dimension ", dim(se)[[1]], " ", dim(se)[[2]])
      se <- ScaleData(se, layer = "counts")
      se <- FindVariableFeatures(se)
      se <- RunPCA(se, features = VariableFeatures(se))
      se <- RunUMAP(se, features = VariableFeatures(se))

      # Save data in obj
      obj@se <- se
      rm(se)
      obj@se.org <- obj@se
      message("LOG: Seurat object created")

      obj <- obj_updateParams(obj,
        updateCurrent = TRUE,
        data_path = data_path,
        TEST = test,
        HVF = HVF,
        test_genes = test_genes,
        test_samples = test_samples
      )
    } else {
      message("Compy from original se")
      obj@se <- obj@se.org
    }

    if (!is.null(pathw)) {
      obj <- obj_updateParams(obj, updateCurrent = TRUE, pathw = pathw)
      obj <- obj_setGenes(obj, pathw)
      message("LOG: Number of genes: ", length(obj@params$genes))
      gene_names <- rownames(obj@se)
      gene.flag <- gene_names %in% obj@params$genes
      message("LOG: intersection pathw and seurat: ", sum(gene.flag))
      obj@se <- obj@se[obj@params$genes, ]
    }


    message("Completed Loading")
    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Melanoma", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})

### obj_geMatrixHVF ----
setMethod("obj_getMatrixHVF", "Melanoma", function(obj) {
  m <- obj@se@assays$RNA@layers$counts[which(obj@se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
  m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
  obj@data$m <- m
  return(obj)
})

# . ##############################################################################
# Exp genera ----
## . ##############################################################################
setClass("Exp",
  contains = "database"
)

### obj_createSeuratObject ----
# function(object, data, gene_names, cell_metadata, where.cell_names, pathw, test = FALSE, HVF = FALSE) {
setMethod("obj_createSeuratObject", "Exp", function(obj, se, gene_names, cell_metadata, where.cell_names, pathw = NULL, test = FALSE, HVF = FALSE, test_genes = 300, test_samples = 500) {
  if (length(where.cell_names) == 2) {
    new.names <- paste0(cell_metadata[[where.cell_names[1]]], "_", cell_metadata[[where.cell_names[2]]])
    cell_metadata$new.names <- new.names
  } else if (length(where.cell_names) == 1) {
    new.names <- cell_metadata[[where.cell_names[1]]]
  } else {
    stop("Invalid where.cell_names")
  }

  if (!is.null(pathw)) {
    obj@params$pathw <- pathw
    obj <- obj_setGenes(obj, pathw)
    gene.flag <- gene_names %in% obj@genes
    se <- se[gene.flag, ]
    gene_names <- gene_names[gene.flag]
  }


  if (test) {
    # tgenes=nrow(se)
    # se <- se[1:tgenes, 1:tsamples]
    tsamples <- min(test_samples, ncol(se))
    se <- se[, 1:tsamples]
    cell_metadata <- cell_metadata[1:tsamples, ]
    gene_names <- gene_names[Matrix::rowSums(se) > 0]
    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
  }
  if (debug) {
    message("Dimension matrix is ", dim(se)[[1]], " ", dim(se)[[2]])
    message("Dimension gene_names ", length(gene_names), " cell_metadata ", dim(cell_metadata)[[1]])
  }


  # newmax=1200*1024^2
  # options(future.globals.maxSize=newmax)

  obj@se <- CreateSeuratObject(counts = se)
  message("LOG: Seurat object has dimension ", dim(se)[[1]], " ", dim(se)[[2]])
  obj@se <- ScaleData(obj@se, layer = "counts")
  obj@se <- FindVariableFeatures(obj@se)
  obj@se <- RunPCA(obj@se, features = VariableFeatures(obj@se))
  obj@se <- RunUMAP(obj@se, features = VariableFeatures(obj@se))

  # dimnames(obj@se) <- list(gene_names, new.names)
  obj@se@assays$RNA@layers$counts@Dimnames <- list(gene_names, cell_metadata$new.names)
  obj@se@meta.data <- cell_metadata

  # if (HVF) {
  #  m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
  # } else {
  #  data <- se@assays$RNA@layers$counts
  # }

  return(obj)
})

# . ##############################################################################
# Exp1 ----
## . ##############################################################################
setClass("Exp1",
  contains = "Exp"
)

### obj_loadData ----
setMethod(
  "obj_loadData", "Exp1",
  function(obj,
           data_path = "/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500, ...) {
    # isnew <- obj_areParamsEqual(obj, update = TRUE)
    if (FALSE) {
      data_path <- "/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx"
    }
    if (is.null(obj@se.org)) { # | isnew) {
      message("LOG: Full loading")
      cell_metadata <- read.table("/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt", header = TRUE, sep = "\t")
      gene_names <- read.table("/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")
      se <- Matrix::readMM(data_path)
      se <- t(se)

      obj <- obj_createSeuratObject(obj, se, gene_names$V1, cell_metadata, where.cell_names = c("Library", "Cell.barcode"), pathw, test, HVF, test_genes, test_samples)
      obj@se.org <- obj@se
      message("LOG: Seurat object created")
    } else {
      message("Compy from original se")
      obj@se <- obj@se.org
    }

    if (!is.null(pathw)) {
      obj <- obj_updateParams(obj, updateCurrent = TRUE, pathw = pathw)
      obj <- obj_setGenes(obj, pathw)
      message("LOG: Number of genes: ", length(obj@params$genes))
      genenames <- rownames(obj@se)
      gene.flag <- genenames %in% obj@params$genes
      message("LOG: intersection pathw and seurat: ", sum(gene.flag))
      obj@se <- obj@se[gene.flag, ]
      gene_names <- gene_names[gene.flag]
    }

    # obj <- obj_updateParams(obj,
    #   updateCurrent = TRUE,
    #   data_path = data_path,
    #   TEST = test,
    #   HVF = HVF
    # )
    message("Completed Loading")
    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Exp1", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})

# . #############################################################################
# Exp2 ----
## . ############################################################################
# Define the 'Exp2' Class that inherits from 'database'
setClass("Exp2",
  contains = "Exp"
)

### obj_loadData ----
setMethod(
  "obj_loadData", "Exp2",
  function(obj,
           data_path = "/app/data/AllonKleinLab/Experiment2/stateFate_inVivo_normed_counts.mtx",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500, ...) {
    # isnew <- obj_areParamsEqual(obj, update = TRUE)
    if (is.null(obj@se.org)) { #| isnew) {
      message("LOG: Full loading")
      gene_names <- read.table("/app/data/AllonKleinLab/Experiment2/stateFate_inVivo_gene_names.txt", sep = "\t")
      cell_metadata <- read.table("/app/data/AllonKleinLab/Experiment2/stateFate_inVivo_metadata.txt", header = TRUE, sep = "\t")

      se <- Matrix::readMM(data_path)
      se <- t(se)

      obj <- obj_createSeuratObject(obj, se, gene_names$V1, cell_metadata, c("Library", "Cell.barcode"), pathw)
      message("LOG: Seurat object created")

      obj@se.org <- obj@se
    } else {
      obj@se <- obj@se.org
    }

    if (!is.null(pathw)) {
      obj <- obj_updateParams(obj, updateCurrent = TRUE, pathw = pathw)
      obj <- obj_setGenes(obj, pathw)
      message("LOG: Number of genes: ", length(obj@params$genes))
      gene.flag <- gene_names %in% obj@genes
      message("LOG: intersection pathw and seurat: ", sum(gene.flag))
      se <- se[gene.flag, ]
      gene_names <- gene_names[gene.flag]
    }

    # obj <- obj_updateParams(obj,
    #   updateCurrent = TRUE,
    #   data_path = data_path,
    #   TEST = test,
    #   HVF = HVF,
    #   test_genes = test_genes,
    #   test_samples = test_samples
    # )
    message("Completed Loading")
    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Exp2", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})

# . #############################################################################
# Exp3 ----
## . ############################################################################
# Define the 'Exp3' Class that inherits from 'database'
setClass("Exp3",
  contains = "Exp"
)

### obj_loadData ----
setMethod(
  "obj_loadData", "Exp3",
  function(obj,
           data_path = "/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_normed_counts.mtx",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500,
           ...) {
    if (FALSE) {
      data_path <- "/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_normed_counts.mtx"
    }

    # isnew <- obj_areParamsEqual(obj, update = TRUE)
    if (is.null(obj@se.org)) { # | isnew) {
      message("LOG: Full loading")
      gene_names <- read.table("/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_gene_names.txt", sep = "\t")
      cell_metadata <- read.table("/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_metadata.txt", header = TRUE, sep = "\t")

      se <- Matrix::readMM(data_path)
      se <- t(se)
      obj <- obj_createSeuratObject(obj, se, gene_names$V1, cell_metadata, where.cell_names = c("Library", "Cell.barcode"), pathw)
      message("LOG: Seurat object created")

      obj@se.org <- obj@se
    } else {
      obj@se <- obj@se.org
    }

    if (!is.null(pathw)) {
      obj <- obj_updateParams(obj, updateCurrent = TRUE, pathw = pathw)
      obj <- obj_setGenes(obj, pathw)
      message("LOG: Number of genes: ", length(obj@params$genes))
      gene.flag <- gene_names %in% obj@genes
      message("LOG: intersection pathw and seurat: ", sum(gene.flag))
      se <- se[gene.flag, ]
      gene_names <- gene_names[gene.flag]
    }

    # obj <- obj_updateParams(obj,
    #  updateCurrent = TRUE,
    #  data_path = data_path,
    #  TEST = test,
    #  HVF = HVF,
    #  test_genes = test_genes,
    #  test_samples = test_samples
    # )

    message("Completed Loading")
    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Exp3", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})

# . #############################################################################
# MouseCortex ----
## . ############################################################################
setClass("MouseCortex",
  contains = "database"
)

### obj_loadData ----
setMethod(
  "obj_loadData", "MouseCortex",
  function(obj,
           data_path = "/app/data/MouseCortex/MouseCortex.RData",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500,
           ...) {
    load(data_path)
    mv("MouseCortex", "se")

    raw_counts <- se@raw.data
    normalized_data <- se@data
    scaled_data <- se@scale.data
    var_genes <- se@var.genes
    meta_data <- se@meta.data
    se.ident <- se@ident
    rm(se)

    se <- CreateSeuratObject(counts = raw_counts, meta.data = meta_data)
    rm(raw_counts, meta_data)
    se[["RNA"]] <- SetAssayData(se[["RNA"]], layer = "data", new.data = normalized_data)
    rm(normalized_data)
    se[["RNA"]] <- SetAssayData(se[["RNA"]], layer = "scale.data", new.data = scaled_data)
    rm(scaled_data)
    VariableFeatures(se) <- var_genes
    rm(var_genes)
    Idents(se) <- se.ident
    rm(se.ident)

    if (!is.null(pathw)) {
      obj@params$pathw <- pathw
      obj <- obj_setGenes(obj, pathw)
      se <- se[obj@params$genes, ]
    } else if (test) {
      tgenes <- min(test_genes, nrow(se))
      tsamples <- min(test_samples, ncol(se))
      se <- se[1:tgenes, 1:tsamples]
      rm(tgenes, tsamples)
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    }

    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    se <- FindVariableFeatures(se)
    se <- RunPCA(se, features = VariableFeatures(se))
    se <- RunTSNE(se)
    se <- RunUMAP(se, features = VariableFeatures(se))

    if (HVF) {
      m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
    } else {
      m <- se@assays$RNA@layers$counts
    }
    m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
    m <- as.matrix(m)

    obj@se <- se
    obj@data$m <- m

    obj@params$data_path <- data_path
    obj@params$TEST <- test
    obj@params$HVF <- HVF
    obj@params$test_genes <- test_genes
    obj@params$test_samples <- test_samples

    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "MouseCortex", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})

# . #############################################################################
# Myocardial ----
## . ############################################################################
setClass("Myocardial",
  contains = "database"
)

### obj_loadData ----
setMethod(
  "obj_loadData", "Myocardial",
  function(obj,
           data_path = "/app/data/MyocardialInfarction/e61af320-303a-4029-8500-db6636bba0d4.rds",
           test = FALSE,
           HVF = TRUE,
           test_genes = 300,
           test_samples = 500,
           ...) {
    se <- readRDS(data_path)

    if (!is.null(pathw)) {
      obj@params$pathw <- pathw
      obj <- obj_setGenes(obj, pathw)
      se <- se[obj@params$genes, ]
    } else if (test) {
      tgenes <- min(test_genes, nrow(se))
      tsamples <- min(test_samples, ncol(se))
      se <- se[1:tgenes, 1:tsamples]
      rm(tgenes, tsamples)
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    }

    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    se <- ScaleData(se, layer = "counts")
    se <- FindVariableFeatures(se)

    se <- RunPCA(se)
    se <- RunUMAP(se, features = VariableFeatures(se))

    if (HVF) {
      m <- se@assays$RNA@counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
    } else {
      m <- se@assays$RNA@layers$counts
    }
    m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
    m <- as.matrix(m)

    obj@se <- se
    obj@data$m <- m

    obj@params$data_path <- data_path
    obj@params$TEST <- test
    obj@params$HVF <- HVF
    obj@params$test_genes <- test_genes
    obj@params$test_samples <- test_samples

    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Myocardial", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})
