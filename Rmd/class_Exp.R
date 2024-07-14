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

  # Creating new names for removing duplicated
  if (length(where.cell_names) == 2) {
    new.names <- paste0(cell_metadata[[where.cell_names[1]]], "_", cell_metadata[[where.cell_names[2]]])
    cell_metadata$new.names <- new.names
  } else if (length(where.cell_names) == 1) {
    new.names <- cell_metadata[[where.cell_names[1]]]
  } else {
    stop("Invalid where.cell_names")
  }

  if (debug) message("DEBUG: data matrix has dimension post rem0 ", dim(se)[[1]], " ", dim(se)[[2]])

  # REMOVE duplicates
  se <- se[!duplicated(gene_names), !duplicated(cell_metadata$new.names)]
  gene_names <- gene_names[!duplicated(gene_names)]
  cell_metadata <- cell_metadata[!duplicated(cell_metadata$new.names), ]

  # CREATE OBJECT
  obj@se <- CreateSeuratObject(counts = se) # , meta.data = cell_metadata)
  if (debug) message("DEBUG: Seurat object has dimension ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])

  rownames(obj@se) <- gene_names
  colnames(obj@se) <- cell_metadata$new.names
  obj@se <- AddMetaData(obj@se, metadata = cell_metadata)

  # SAVE obj@se to obj@se.org
  obj@se.org <- obj@se

  # TEST
  if (obj@params$test) {
    if (!is.null(obj@params$pathw)) {
      tgenes <- nrow(obj@se)
    } else {
      tgenes <- min(obj@params$test_genes, nrow(obj@se))
    }

    tsamples <- min(obj@params$test_samples, ncol(obj@se))

    obj@se <- obj@se[1:tgenes, 1:tsamples]
    # row_filter <- Matrix::rowSums(obj@se) > 0
    # col_filter <- Matrix::colSums(obj@se) > 0
    # obj@se <- obj@se[row_filter, col_filter]
    obj@se <- obj@se[Matrix::rowSums(obj) > 0, Matrix::colSums(obj_getSeData(obj)) > 0]
    if (debug) message("DEBUG: Seurat after test has dimension ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])
  }
  # obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
  # obj@se <- FindVariableFeatures(obj@se)
  # obj@se <- RunPCA(obj@se, features = rownames(obj@se))
  # obj@se <- RunUMAP(obj@se, features = rownames(obj@se))

  obj@se <- ScaleData(obj@se, layer = "counts")
  obj@se <- FindVariableFeatures(obj@se)
  # obj@se <- RunPCA(obj@se, features = VariableFeatures(obj@se))
  # obj@se <- RunUMAP(obj@se, features = VariableFeatures(obj@se))

  if (obj@params$hvf) {
    if (debug) message("DEBUG: obj_createSeuratObject | number of genes with rank ", sum(which(obj@se@assays$RNA@meta.data$vf_vst_counts_rank > 0)))
    obj@se <- obj@se[which(obj@se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]

    # obj@se <- obj@se[Matrix::rowSums(obj@se) > 0, Matrix::colSums(obj@se) > 0]
    obj@se <- obj@se[Matrix::rowSums(obj_getSeData(obj@se)) > 0, Matrix::colSums(obj@se) > 0]
    message("LOG: HVF: new dimension of se is ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])
    obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
    obj@se <- FindVariableFeatures(obj@se)
    obj@se <- RunPCA(obj@se, features = rownames(obj@se))
    obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
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
    if (!is.null(obj@se.org)) {
      message("LOG: obj_loadData | Copy from original se")
      obj@se <- obj@se.org
    } else {
      message("LOG: Full loading")
      cell_metadata <- read.table("/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt", header = TRUE, sep = "\t")
      gene_names <- read.table("/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")

      se <- Matrix::readMM(data_path)
      se <- t(se)

      obj <- obj_createSeuratObject(obj, se, gene_names = gene_names$V1, cell_metadata = cell_metadata, where.cell_names = c("Library", "Cell.barcode"))

      message("LOG: Seurat object created")
      obj@se.org <- obj@se
    }

    if (!is.null(obj@params$pathw)) {
      # obj <- obj_updateParams(obj, updateCurrent = TRUE, pathw = pathw)
      obj <- obj_setGenes(obj)
      gene.flag <- rownames(obj@se) %in% obj@params$genes
      message("LOG: Number of genes: ", length(obj@params$genes))
      message("LOG: intersection pathw and genenames: ", sum(gene.flag))
      obj@se <- obj@se[gene.flag, ]

      obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
      obj@se <- FindVariableFeatures(obj@se)
      obj@se <- RunPCA(obj@se, features = rownames(obj@se))
      obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
    }

    message("Completed Loading")
    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Exp1", function(obj) {
  return(obj@se@assays$RNA@layers$counts)
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
    if (!is.null(obj@se.org)) { #| isnew) {
      message("LOG: obj_loadData | Copy from original se")
      obj@se <- obj@se.org
    } else {
      message("LOG: Full loading")
      gene_names <- read.table("/app/data/AllonKleinLab/Experiment2/stateFate_inVivo_gene_names.txt", sep = "\t")
      cell_metadata <- read.table("/app/data/AllonKleinLab/Experiment2/stateFate_inVivo_metadata.txt", header = TRUE, sep = "\t")

      se <- Matrix::readMM(data_path)
      se <- t(se)

      obj <- obj_createSeuratObject(obj, se, gene_names$V1, cell_metadata, c("Library", "Cell.barcode"))

      message("LOG: Seurat object created")
      obj@se.org <- obj@se
    }

    if (!is.null(obj@params$pathw)) {
      obj <- obj_setGenes(obj)
      gene.flag <- rownames(obj@se) %in% obj@params$genes
      message("LOG: Number of genes: ", length(obj@params$genes))
      message("LOG: intersection pathw and genenames: ", sum(gene.flag))
      obj@se <- obj@se[gene.flag, ]

      obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
      obj@se <- FindVariableFeatures(obj@se)
      obj@se <- RunPCA(obj@se, features = rownames(obj@se))
      obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
    }

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
# . ############################################################################
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
    if (!is.null(obj@se.org)) {
      obj@se <- obj@se.org
    } else {
      message("LOG: Full loading")
      gene_names <- read.table("/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_gene_names.txt", sep = "\t")
      cell_metadata <- read.table("/app/data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_metadata.txt", header = TRUE, sep = "\t")

      se <- Matrix::readMM(data_path)
      se <- t(se)

      obj <- obj_createSeuratObject(obj, se, gene_names$V1, cell_metadata, where.cell_names = c("Library", "Cell.barcode"))

      message("LOG: Seurat object created")
      obj@se.org <- obj@se
    }

    if (!is.null(obj@params$pathw)) {
      obj <- obj_setGenes(obj)
      gene.flag <- rownames(obj@se) %in% obj@params$genes
      message("LOG: Number of genes: ", length(obj@params$genes))
      message("LOG: intersection pathw and genenames: ", sum(gene.flag))
      obj@se <- obj@se[gene.flag, ]

      obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
      obj@se <- FindVariableFeatures(obj@se)
      obj@se <- RunPCA(obj@se, features = rownames(obj@se))
      obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
    }

    message("Completed Loading")
    return(obj)
  }
)

### obj_getSeData ----
setMethod("obj_getSeData", "Exp3", function(obj) {
  return(obj@se@assays$RNA@layers$counts)
})
