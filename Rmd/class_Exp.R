# . ##############################################################################
# Exp genera ----
## . ##############################################################################
setClass("Exp",
  contains = "database"
)

### obj_createSeuratObject ----
# function(object, data, gene_names, cell_metadata, where.cell_names, pathw, test = FALSE, HVF = FALSE) {
setMethod("obj_createSeuratObject", "Exp", function(obj, se, gene_names, cell_metadata, where.cell_names) {
  if (debug) message("DEBUG: classExp | Entering in function")
  if (length(where.cell_names) == 2) {
    new.names <- paste0(cell_metadata[[where.cell_names[1]]], "_", cell_metadata[[where.cell_names[2]]])
    cell_metadata$new.names <- new.names
  } else if (length(where.cell_names) == 1) {
    new.names <- cell_metadata[[where.cell_names[1]]]
  } else {
    stop("Invalid where.cell_names")
  }

  if (debug) message("DEBUG: data matrix has dimension post rem0 ", dim(se)[[1]], " ", dim(se)[[2]])

  if (obj@params$test) {
    if (!is.null(obj@params$pathw)) {
      tgenes <- nrow(se)
    } else {
      tgenes <- min(obj@params$test_genes, nrow(se))
    }

    tsamples <- min(obj@params$test_samples, ncol(se))
    se <- se[1:tgenes, 1:tsamples]
    gene_names <- gene_names[1:tgenes]
    cell_metadata <- cell_metadata[1:tsamples, ]

    row_filter <- Matrix::rowSums(se) > 0
    col_filter <- Matrix::colSums(se) > 0

    se <- se[row_filter, col_filter]
    gene_names <- gene_names[row_filter]
    cell_metadata <- cell_metadata[col_filter, ]
  }

  obj@se <- CreateSeuratObject(counts = se)
  if (debug) message("DEBUG: Seurat object has dimension ", dim(se)[[1]], " ", dim(se)[[2]])
  obj@se <- ScaleData(obj@se, layer = "counts")
  obj@se <- FindVariableFeatures(obj@se)
  obj@se <- RunPCA(obj@se, features = VariableFeatures(obj@se))
  obj@se <- RunUMAP(obj@se, features = VariableFeatures(obj@se))

  # dimnames(obj@se) <- list(gene_names, new.names)
  obj@se@assays$RNA@layers$counts@Dimnames <- list(gene_names, cell_metadata$new.names)
  obj@se@meta.data <- cell_metadata

  if (obj@params$hvf) {
    obj@se <- obj@se[which(obj@se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
    obj@se <- obj@se[Matrix::rowSums(obj@se) > 0, Matrix::colSums(obj@se) > 0]
    message("LOG: HVF: new dimension of se is ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])
  }

  # if (!is.null(pathw)) {
  #  obj@params$pathw <- pathw
  #  obj <- obj_setGenes(obj, pathw)
  #  message("LOG: Number of genes: ", length(obj@params$genes))
  #  gene.flag <- gene_names %in% obj@genes
  #  message("LOG: intersection pathw and seurat: ", sum(gene.flag))
  #  se <- se[gene.flag, ]
  #  gene_names <- gene_names[gene.flag]
  # }

  obj@se.org <- obj@se
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
           ...) {
    # isnew <- obj_areParamsEqual(obj, update = TRUE)
    if (FALSE) {
      data_path <- "/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx"
    }
    if (is.null(obj@se.org)) {
      message("LOG: Full loading")
      cell_metadata <- read.table("/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt", header = TRUE, sep = "\t")
      gene_names <- read.table("/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")
      se <- Matrix::readMM(data_path)
      se <- t(se)

      obj <- obj_createSeuratObject(obj, se, gene_names$V1, cell_metadata, where.cell_names = c("Library", "Cell.barcode"), obj@params$pathw, obj@params$test, obj@params$hvf, obj@params$test_genes, obj@params$test_samples)
      obj@se.org <- obj@se
      message("LOG: Seurat object created")
    } else {
      message("Compy from original se")
      obj@se <- obj@se.org
    }

    if (!is.null(obj@params$pathw)) {
      # obj <- obj_updateParams(obj, updateCurrent = TRUE, pathw = pathw)
      obj <- obj_setGenes(obj)
      message("LOG: Number of genes: ", length(obj@params$genes))
      genenames <- rownames(obj@se)
      gene.flag <- genenames %in% obj@params$genes
      message("LOG: intersection pathw and seurat: ", sum(gene.flag))
      obj@se <- obj@se[gene.flag, ]
      # gene_names <- gene_names[gene.flag]
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
           ...) {
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

    if (!is.null(obj@params$pathw)) {
      obj <- obj_setGenes(obj)
      message("LOG: Number of genes: ", length(obj@params$genes))
      gene.flag <- gene_names %in% obj@genes
      message("LOG: intersection pathw and seurat: ", sum(gene.flag))
      se <- se[gene.flag, ]
      # gene_names <- gene_names[gene.flag]
    }

    obj <- obj_updateParams(obj,
      updateCurrent = TRUE,
      data_path = data_path
    )
    message("Completed Loading")
    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Exp2", function(obj) {
  return(obj@se@assays$RNA@layers$counts)
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
           ...) {
    if (FALSE) {
      data_path <- "/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_normed_counts.mtx"
    }

    # isnew <- obj_areParamsEqual(obj, update = TRUE)
    if (is.null(obj@se.org)) {
      message("LOG: Full loading")
      gene_names <- read.table("/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_gene_names.txt", sep = "\t")
      cell_metadata <- read.table("/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_metadata.txt", header = TRUE, sep = "\t")

      se <- Matrix::readMM(data_path)
      se <- t(se)
      obj <- obj_createSeuratObject(obj, se, gene_names$V1, cell_metadata, where.cell_names = c("Library", "Cell.barcode"), obj@params$pathw)
      message("LOG: Seurat object created")

      obj@se.org <- obj@se
    } else {
      obj@se <- obj@se.org
    }

    if (!is.null(obj@params$pathw)) {
      obj <- obj_setGenes(obj)
      message("LOG: Number of genes: ", length(obj@params$genes))
      gene.flag <- gene_names %in% obj@genes
      message("LOG: intersection pathw and seurat: ", sum(gene.flag))
      se <- se[gene.flag, ]
      # gene_names <- gene_names[gene.flag]
    }

    obj <- obj_updateParams(obj,
      updateCurrent = TRUE,
      data_path = data_path
    )

    message("Completed Loading")
    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Exp3", function(obj) {
  return(obj@se@assays$RNA@layers$counts)
})
