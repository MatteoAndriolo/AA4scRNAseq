# GENERIC --------
# Define the 'database' Class
setClass("database",
  slots = list(
    se = "Seurat",
    a = "ANY", # Archetype model will be stored here
    m = "matrix", # Data matrix
    umap.archetypes = "ANY",
    combined_plot = "ANY",
    elbowplot = "ANY",
    umap = "ANY",
    pca = "ANY",
    simplexplot = "ANY",
    pathw = "ANY",
    genes = "ANY"
  )
)


## obj_loadData ----
# Generic method for loading data
setGeneric("obj_loadData", function(obj,
                                    data_path = NULL,
                                    test = FALSE,
                                    HVF = TRUE,
                                    pathw = NULL,
                                    test_genes = 300,
                                    test_samples = 500) {
  standardGeneric("obj_loadData")
})

## obj_getSeData ------------------------------------------------------------
setGeneric("obj_getSeData", function(obj) {
  standardGeneric("obj_getSeData")
})

## obj_visualizeData -----------------------------------------------------------
# Method to visualize data
setGeneric("obj_visualizeData", function(obj, out_path) {
  standardGeneric("obj_visualizeData")
})

setMethod("obj_visualizeData", "database", function(obj, out_path) {
  se <- obj@se

  imgname_pca <- sprintf("%s/PCA.png", out_path)
  message(sprintf("Saving Image --- %s", imgname_pca))
  pcaplot <- PCAPlot(se)
  obj@pca <- pcaplot

  imgname_umap <- sprintf("%s/UMAP.png", out_path)
  message(sprintf("Saving Image --- %s", imgname_umap))
  umapplot <- UMAPPlot(se)
  umapplot
  obj@umap <- umapplot

  # Extract legends
  pca_legend <- cowplot::get_legend(pcaplot)
  umap_legend <- cowplot::get_legend(umapplot)

  # Combine plots without legends
  combined_plot <- plot_grid(
    pcaplot + theme(legend.position = "none"),
    umapplot + theme(legend.position = "none"),
    labels = c("A", "B"),
    ncol = 2
  )

  # Combine legends
  # combined_legend <- plot_grid(pca_legend, umap_legend, ncol = 1)
  # Onli one legend
  combined_legend <- plot_grid(umap_legend, ncol = 1)

  # Combine plots and legends
  final_plot <- plot_grid(
    combined_plot, combined_legend,
    ncol = 2, rel_widths = c(6, 1)
  )

  obj@combined_plot <- final_plot
  print(final_plot)


  # combined_plot <-
  #   plot_grid(pcaplot, umapplot, labels = c("A", "B"))
  # obj@combined_plot <- combined_plot
  # print(combined_plot)


  elbowplot <- ElbowPlot(se)
  elbowplot
  obj@elbowplot <- elbowplot
  print(elbowplot)

  if (!is.null(out_path)) {
    imgname <- sprintf("%s/Combined_Plots.png", out_path)
    message(sprintf("Saving Image --- %s", imgname))
    ggsave(imgname, plot = combined_plot)

    imgname <- sprintf("%s/elbow.pdf", out_path)
    message(sprintf("Saving Imagine --- %s", imgname))
    ggsave(imgname, plot = elbowplot)
  }
  return(obj)
})

## obj_performArchetypes -------
# Method to perform archetypal analysis
setGeneric("obj_performArchetypes", function(obj, k = 5, HVF = TRUE) {
  standardGeneric("obj_performArchetypes")
})
setMethod("obj_performArchetypes", "database", function(obj, k = 5, HVF = TRUE) {
  m <- obj@m
  m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
  m <- as.matrix(m)

  # Time this and print as message
  tstart <- Sys.time()
  a <- tryCatch(
    {
      archetypes::archetypes(
        m,
        k = k,
        verbose = TRUE,
        maxIterations = 10,
        saveHistory = TRUE
      )
    },
    error = function(e) {
      message(sprintf("Error in archetypes computation: %s", e$message))
      return(NULL)
    }
  )
  tend <- Sys.time()
  message(sprintf("Archetypes Computed in %s", tend - tstart))

  if (!is.null(a)) {
    obj@a <- a
  } else {
    stop("Archetypes computation failed, obj@a not assigned.")
  }

  save(a, file = sprintf("%s/Archetypes_%02d.rds", out_path, k))
  return(obj)
})
# setMethod("obj_performArchetypes", "database", function(obj, k = 5, HVF = TRUE) {
#  m <- obj@m
#  m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
#  m <- as.matrix(m)
#
#  #time this and print as message
#  tstart <- Sys.time()
#  a <-
#    archetypes::archetypes(
#      m,
#      k = k,
#      verbose = TRUE,
#      maxIterations = 10,
#      saveHistory = TRUE
#    )
#  tend <- Sys.time()
#  message(sprintf("Archetypes Computed in %s", tend - tstart))
#  obj@a <- a
#
#  return(obj)
# })

## obj_visualizeArchetypes -----------------------------------------------------
# Method to visualize archetypes
setGeneric("obj_visualizeArchetypes", function(obj, out_path) {
  standardGeneric("obj_visualizeArchetypes")
})

setMethod("obj_visualizeArchetypes", "database", function(obj, out_path) {
  a <- obj@a
  k <- a$k
  m <- obj@m
  imgname <- sprintf("%s/Archetypes_%2d.png", out_path, k)
  plotarchetyps <- xyplot(a, as.matrix(obj_getSeData(obj)))
  plotarchetyps
  obj@simplexplot <- plotarchetyps
  return(obj)
})

## obj_umapArchetypes ----------------------------------------------------------
setGeneric("obj_umapArchetypes", function(obj, out_path = NULL, treshold = 0.2) {
  standardGeneric("obj_umapArchetypes")
})
setMethod("obj_umapArchetypes", "database", function(obj,
                                                     out_path = NULL,
                                                     treshold = 0.2) {
  se <- obj@se
  a <- obj@a
  k <- a$k

  umap_result <- UMAPPlot(se)
  umap_data <- as.data.frame(umap_result$data)[, 1:2]
  colnames(umap_data) <- c("UMAP1", "UMAP2")

  # Get the archetype weights
  weights <- coef(a)
  weights <- as.data.frame(weights)
  # Set a minimum threshold
  weights[weights < treshold] <- 0

  column_sums <- colSums(a$archetypes) # Sum of each column
  normalized_mat <-
    sweep(a$archetypes, 2, column_sums, FUN = "/") # Divide each element by its column sum
  weights <- as.data.frame(normalized_mat)

  plot_list <- list()
  # Plotting
  for (i in 1:k) {
    umap_data$weight <- t(weights[i, ])
    plot_title <- sprintf("UMAP Archetype %d", i)

    umap_plot <-
      ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = weight)) +
      geom_point(size = 1) +
      scale_color_gradient(low = "grey", high = "red") +
      ggtitle(plot_title) +
      labs(color = "Weight")
    print(umap_plot)

    if (!is.null(out_path)) {
      imgname <- sprintf("%s/UMAP_Archetype_%d.%d.png", out_path, k, i)
      ggsave(imgname, plot = umap_plot)
      message(sprintf("Saving Image --- %s", imgname))
    }

    plot_list[[i]] <- umap_plot
  }

  # Combine all plots into a single image
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)

  # Save the combined image
  print(combined_plot)
  if (!is.null(out_path)) {
    combined_imgname <- sprintf("%s/UMAP_Combined.png", out_path)
    ggsave(
      combined_imgname,
      plot = combined_plot,
      width = 12,
      height = 8
    )
    message(sprintf("Saving Combined Image --- %s", combined_imgname))
  }
  obj@umap.archetypes <- combined_plot
  return(obj)
})

## obj_setGenes ---------------------------------------------------------------
setGeneric("obj_setGenes", function(obj, pathw, pathGenes = "/app/data/list_genes_pathway.RData") {
  standardGeneric("obj_setGenes")
})

setMethod("obj_setGenes", "database", function(obj, pathw, pathGenes = "/app/data/list_genes_pathway.RData") {
  load(pathGenes)
  list_genes_human_pathway <- list_genes_human_pathway
  list_genes_mouse_pathway <- list_genes_mouse_pathway

  message("Setting Genes %s", pathw)
  if (inherits(obj, "Melanoma")) {
    if (is.character(pathw) && length(pathw) == 1) {
      genes <- list_genes_human_pathway[[pathw]]
    } else if (is.list(pathw) || is.vector(pathw)) {
      genes <- lapply(pathw, function(p) list_genes_human_pathway[[p]])
    }
  } else {
    if (is.character(pathw) && length(pathw) == 1) {
      genes <- list_genes_mouse_pathway[[pathw]]
    } else if (is.list(pathw) || is.vector(pathw)) {
      genes <- lapply(pathw, function(p) list_genes_mouse_pathway[[p]])
    }
  }

  if (is.null(genes)) {
    stop("Pathway not found")
  }

  obj@pathw <- pathw
  obj@genes <- genes
  return(obj)
})

# ##############################################################################
# Melanoma ---------------------------------------------------------------------
# ##############################################################################
setClass("Melanoma",
  contains = "database"
)

## obj_loadData --------------------------------------------------------------------
setMethod(
  "obj_loadData", "Melanoma",
  function(obj,
           data_path = "../data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500) {
    se <- read.table(data_path, header = TRUE)
    se <- se[!duplicated(se[, 1]), ] # remove duplicated genes
    rownames(se) <- se[, 1] # extract from matrix rownames
    se <-
      se[, 2:ncol(se)] # elide rownames from gene expression matrix

    metadata <- se[1:3, ] # extract metadata
    metadata <- t(metadata) %>%
      data.frame() %>%
      mutate(across(where(is.character), as.numeric))

    se <- se[4:nrow(se), ]
    se <- se %>%
      data.frame() %>%
      mutate(across(where(is.character), as.numeric))
    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]


    if (test) {
      message("No Pathway")
      tgenes <- min(test_genes, nrow(se))
      tsamples <- min(test_samples, ncol(se))
      metadata <- metadata[1:tsamples, ]
      se <- se[1:tgenes, 1:tsamples]
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    }

    # se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]

    se <- CreateSeuratobj(counts = se, meta.data = metadata)
    se <- ScaleData(se, layer = "counts")
    se <- FindVariableFeatures(se)
    se <- RunPCA(se, features = VariableFeatures(se))
    se <- RunUMAP(se, features = VariableFeatures(se))
    obj@se <- se

    if (!is.null(pathw)) {
      message("Pathway")
      obj <- obj_setGenes(obj, pathw)
      message(length(obj@genes))
      se <- se[obj@genes, ]
    }

    if (HVF) {
      m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
    } else {
      m <- se@assays$RNA@layers$counts
      m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
    }
    message("finished HVF")

    m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
    m <- as.matrix(m)
    obj@m <- m

    message("Completed Loading")
    return(obj)
  }
)

## obj_getSeData ------------------------------------------------------------
setMethod("obj_getSeData", "Melanoma", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})

# melanoma = new("Melanoma")
# obj_loadData(
#  melanoma,
#  data_path = here(
#    "/app/data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt"
#  ),
#  test = FALSE,
#  HVF = TRUE,
#  test_genes = 300,
#  test_samples = 500
# )
# obj_visualizeData(melanoma, out_path = here("/app/out/Melanoma/"))
# obj_performArchetypes(melanoma, k = 5, HVF= TRUE)

# ##############################################################################
# Exp1 -------------------------------------------------------------------------
# ##############################################################################
setClass("Exp1",
  contains = "database"
)

## obj_loadData --------------------------------------------------------------------
setMethod(
  "obj_loadData", "Exp1",
  function(obj,
           data_path = "../data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500) {
    # _ # Binary matrix indicating clonal membership of each cell
    # _ # The rows of this file represent cells and correspond to the rows of _counts_matrix_in_vitro_ (above).
    # _ # The columns represent clones. Not every cell belongs to a clone.
    clone_matrix <- Matrix::readMM("../data/AllonKleinLab/Experiment1/stateFate_inVitro_clone_matrix.mtx")
    # List of cells belonging to the neutrophil/monocyte trajectory that were used in becnmark analysis
    neutrophil_monocyte_trajectory <- read.table("../data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_monocyte_trajectory.txt", header = TRUE, sep = "\t")
    # pseudotime for neutrophil trajectory cells
    neutrophil_pseudotime <- read.table("../data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_pseudotime.txt", header = TRUE, sep = "\t")
    cell_metadata <- read.table("../data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt", header = TRUE, sep = "\t")
    gene_names <- read.table("../data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")


    se <- Matrix::readMM(data_path)
    se@Dimnames[[2]] <- gene_names$V1
    # se@Dimnames[[1]]=cell_metadata$Cell.barcode
    # remove duplicated rows columns looking at names
    se <- se[!duplicated(se@Dimnames[[1]]), !duplicated(se@Dimnames[[2]])]
    # se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    se <- CreateSeuratObject(counts = se)
    se@assays$RNA@layers$counts@Dimnames[[2]] <- gene_names$V1

    if (!is.null(pathw)) {
      obj@pathw <- pathw
      obj <- obj_setGenes(obj, pathw)
      se <- se[gene_names$V1, ]
    } else if (test) {
      tgenes <- min(TEST_genes, nrow(se))
      tsamples <- min(TEST_samples, ncol(se))
      se <- se[1:tgenes, 1:tsamples]
      rm(tgenes, tsamples)
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    }

    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    se <- ScaleData(se, layer = "counts")
    se <- FindVariableFeatures(se)

    se <- RunPCA(se, features = VariableFeatures(se))
    se <- RunUMAP(se, features = VariableFeatures(se))

    if (HVF) {
      m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
    } else {
      m <- se@assays$RNA@layers$counts
    }

    m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
    m <- as.matrix(m)

    obj@se <- se
    obj@m <- m

    return(obj)
  }
)

## obj_getSeData ------------------------------------------------------------
setMethod("obj_getSeData", "Exp1", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})


# ##############################################################################
# Exp2 -------------------------------------------------------------------------
# ##############################################################################
# Define the 'Melanoma' Class that inherits from 'database'
setClass("Exp2",
  contains = "database"
)

## obj_loadData --------------------------------------------------------------------
setMethod(
  "obj_loadData", "Exp2",
  function(obj,
           data_path = "../data/AllonKleinLab/Experiment2/stateFate_inVivo_normed_counts.mtx",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500) {
    # TODO implement
    # _ clone_matrix <- Matrix::readMM("data/AllonKleinLab/Experiment2/stateFate_inVivo_clone_matrix.mtx")
    # _ gene_names <- read.table("data/AllonKleinLab/Experiment2/stateFate_inVivo_gene_names.txt",sep="\t")
    # _ metadata <- read.table("data/AllonKleinLab/Experiment2/stateFate_inVivo_metadata.txt",sep="\t")

    # out_path <- "../out/AllonKleinLab/Experiment2"

    se <- Matrix::readMM(data_path)
    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]

    if (!is.null(pathw)) {
      obj@pathw <- pathw
      obj <- obj_setGenes(obj, pathw)
      se <- se[obj@genes, ]
    } else if (test) {
      tgenes <- min(TEST_genes, nrow(se))
      tsamples <- min(TEST_samples, ncol(se))
      se <- se[1:tgenes, 1:tsamples]
      rm(tgenes, tsamples)
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    }

    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    se <- CreateSeuratobj(counts = se)
    se <- ScaleData(se, layer = "counts")
    se <- FindVariableFeatures(se)

    se <- RunPCA(se, features = VariableFeatures(se))
    se <- RunUMAP(se, features = VariableFeatures(se))


    if (HVF) {
      m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
    } else {
      m <- se@assays$RNA@layers$counts
    }

    m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
    m <- as.matrix(m)

    obj@se <- se
    obj@m <- m

    return(obj)
  }
)

## obj_getSeData ------------------------------------------------------------
setMethod("obj_getSeData", "Exp2", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})


# ##############################################################################
# Exp3 -------------------------------------------------------------------------
# ##############################################################################
# Define the 'Melanoma' Class that inherits from 'database'
setClass("Exp3",
  contains = "database"
)

## obj_loadData --------------------------------------------------------------------
setMethod(
  "obj_loadData", "Exp3",
  function(obj,
           data_path = "../data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_normed_counts.mtx",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500) {
    # TODO implement
    # _ data_clone_matrix <- Matrix::readMM("../data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_clone_matrix.mtx")
    data_gene_names <- read.table("../data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_gene_names.txt", sep = "\t")
    data_gene_names <- lapply(data_gene_names, toupper)
    data_metadata <- read.table("../data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_metadata.txt", sep = "\t")

    # out_path <- "../out/AllonKleinLab/Experiment3"

    se <- Matrix::readMM(data_path)
    se@Dimnames[[1]] <- data_gene_names$V1
    se@Dimnames[[2]] <- data_metadata$V2[-1]
    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]

    if (!is.null(pathw)) {
      obj@pathw <- pathw
      obj <- obj_setGenes(obj, pathw)
      common_genes <- intersect(obj@genes, rownames(se))
      # se <- se[obj@genes, ]
      se <- se[common_genes, ]
    } else if (test) {
      tgenes <- min(TEST_genes, nrow(se))
      tsamples <- min(TEST_samples, ncol(se))
      se <- se[1:tgenes, 1:tsamples]
      rm(tgenes, tsamples)
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    }

    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    se <- CreateSeuratobj(counts = se)
    se <- ScaleData(se, layer = "counts")
    se <- FindVariableFeatures(se)

    se <- RunPCA(se, features = VariableFeatures(se))
    se <- RunUMAP(se, features = VariableFeatures(se))

    if (HVF) {
      m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
    } else {
      m <- se@assays$RNA@layers$counts
    }
    m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
    m <- as.matrix(m)

    obj@se <- se
    obj@m <- m

    return(obj)
  }
)

## obj_getSeData ------------------------------------------------------------
setMethod("obj_getSeData", "Exp3", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})

# ##############################################################################
# MouseCortex-------------------------------------------------------------------
# ##############################################################################
# Define the 'Melanoma' Class that inherits from 'database'
setClass("MouseCortex",
  contains = "database"
)

## obj_loadData --------------------------------------------------------------------
setMethod(
  "obj_loadData", "MouseCortex",
  function(obj,
           data_path = "../MouseCortex/MouseCortex.RData",
           test = FALSE,
           HVF = TRUE,
           pathw = NULL,
           test_genes = 300,
           test_samples = 500) {
    # Data definitions and loading
    # out_path <- "../out/MouseCortex"
    load(data_path)
    MouseCortex <- MouseCortex
    mv("MouseCortex", "se")

    # Manually generate Seurat obj S5
    raw_counts <- se@raw.data
    normalized_data <- se@data
    scaled_data <- se@scale.data
    var_genes <- se@var.genes
    meta_data <- se@meta.data
    se.ident <- se@ident
    rm(se)

    se <- CreateSeuratobj(counts = raw_counts, meta.data = meta_data)
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
      obj@pathw <- pathw
      obj <- obj_setGenes(obj, pathw)
      se <- se[obj@genes, ]
    } else if (test) {
      tgenes <- min(TEST_genes, nrow(se))
      tsamples <- min(TEST_samples, ncol(se))
      se <- se[1:tgenes, 1:tsamples]
      rm(tgenes, tsamples)
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    }

    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    se <- FindVariableFeatures(se)
    # se <- RunPCA(se)
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
    obj@m <- m

    return(obj)
  }
)

## obj_getSeData --------------------------------------------------------------------
setMethod("obj_getSeData", "MouseCortex", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})

# ##############################################################################
# Myocardial -------------------------------------------------------------------
# ##############################################################################
setClass("Myocardial",
  contains = "database"
)

## obj_loadData --------------------------------------------------------------------
setMethod(
  "obj_loadData", "Myocardial",
  function(obj,
           data_path = "../data/MyocardialInfarction/e61af320-303a-4029-8500-db6636bba0d4.rds",
           test = FALSE,
           HVF = TRUE,
           test_genes = 300,
           test_samples = 500) {
    # out_path <- "../out/MyocardialInfarction"
    se <- readRDS(data_path)

    if (!is.null(pathw)) {
      obj@pathw <- pathw
      obj <- obj_setGenes(obj, pathw)
      se <- se[obj@genes, ]
    } else if (test) {
      tgenes <- min(TEST_genes, nrow(se))
      tsamples <- min(TEST_samples, ncol(se))
      se <- se[1:tgenes, 1:tsamples]
      rm(tgenes, tsamples)
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    }

    se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
    se <- ScaleData(se, layer = "counts")
    se <- FindVariableFeatures(se)

    se <- RunPCA(se) # , features = VariableFeatures(se))
    se <- RunUMAP(se, features = VariableFeatures(se))

    # TODO fix this
    if (HVF) {
      m <- se@assays$RNA@counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
    } else {
      m <- se@assays$RNA@layers$counts
    }
    m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
    m <- as.matrix(m)

    obj@se <- se
    obj@m <- m

    return(obj)
  }
)

## obj_getSeData --------------------------------------------------------------------
setMethod("obj_getSeData", "Myocardial", function(obj) {
  se <- obj@se
  return(se@assays$RNA@layers$counts)
})
