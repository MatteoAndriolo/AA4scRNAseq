# Exp general ----
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

  if (debug) message("DEBUG: data matrix has dimension", dim(se)[[1]], " ", dim(se)[[2]])

  # REMOVE duplicates
  se <- se[!duplicated(gene_names), !duplicated(cell_metadata$new.names)]
  gene_names <- gene_names[!duplicated(gene_names)]
  cell_metadata <- cell_metadata[!duplicated(cell_metadata$new.names), ]


  # CREATE OBJECT
  obj@se <- CreateSeuratObject(counts = se, meta.data = cell_metadata)
  if (debug) message("DEBUG: Seurat object has dimension ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])

  rownames(obj@se) <- gene_names
  colnames(obj@se) <- cell_metadata$new.names
  obj@se <- AddMetaData(obj@se, metadata = cell_metadata)
  obj@se$ctype <- factor(cell_metadata$Cell.type.annotation, levels = unique(cell_metadata$Cell.type.annotation))
  obj@se <- obj@se[Matrix::rowSums(GetAssayData(obj@se)) > 0, Matrix::colSums(GetAssayData(obj@se)) > 0]
  obj@se <- SetAssayData(obj@se, layer = "scale.data", new.data = as.matrix(GetAssayData(obj@se)))
  if (debug) message("DEBUG: Seurat object has dimension ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])

  if (debug) message("DEBUG: Seurat object has dimension post rem0 post ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])

  # TEST
  if (obj@params$test) {
    if (!is.null(obj@params$pathw)) {
      tgenes <- nrow(obj@se)
    } else {
      tgenes <- min(obj@params$test_genes, nrow(obj@se))
    }

    tsamples <- min(obj@params$test_samples, ncol(obj@se))

    obj@se <- obj@se[1:tgenes, 1:tsamples]
    obj@se <- obj@se[Matrix::rowSums(GetAssayData(obj@se)) > 0, Matrix::colSums(GetAssayData(obj@se)) > 0]
    if (debug) message("DEBUG: Seurat after test has dimension ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])
  }

  if (obj@params$hvf) {
    if (debug) message("DEBUG: obj_createSeuratObject | number of genes with rank ", sum(which(obj@se@assays$RNA@meta.data$vf_vst_counts_rank > 0)))
    obj@se <- obj@se[which(obj@se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]

    # obj@se <- obj@se[Matrix::rowSums(obj@se) > 0, Matrix::colSums(obj@se) > 0]
    obj@se <- obj@se[Matrix::rowSums(GetAssayData(obj@se)) > 0, Matrix::colSums(GetAssayData(obj@se)) > 0]
    message("LOG: HVF: new dimension of se is ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])
    obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
  }

  obj@se <- FindVariableFeatures(obj@se)
  obj@se <- RunPCA(obj@se, features = VariableFeatures(obj@se), seed.use = obj@params$rseed)
  obj@se <- RunUMAP(obj@se, features = VariableFeatures(obj@se), seed.use = obj@params$rseed)
  obj@se <- RunTSNE(obj@se, features = VariableFeatures(obj@se), seed.use = obj@params$rseed)

  str(obj@se)
  obj@se.org <- obj@se
  str(obj@se)
  if (!is.null(obj@params$pathw)) {
    message("LOG: obj_createSeu | setting up pathw ", obj@params$pathw)
    obj <- obj_setGenes(obj)
    gene_names <- rownames(obj@se)
    gene.flag <- gene_names %in% obj@params$genes
    if (debug) {
      message("DEBUG: obj_createSeu | pathw | number genes pathw = ", length(obj@params$genes))
      message("DEBUG: obj_createSeu | pathw | first 5 ", toString(obj@params$genes[1:5]))
      message("DEBUG: obj_createSeu | pathw | number genes data = ", length(gene_names))
      message("DEBUG: obj_createSeu | pathw | first 5 ", toString(gene_names[1:5]))
      message("DEBUG: obj_createSeu | pathw | intersectoin gives ", sum(gene.flag), " genes")
    }
    message("LOG: Number of genes: ", length(obj@params$genes))
    message("LOG: intersection pathw and genenames: ", sum(gene.flag))
    obj@se <- obj@se[gene.flag, ]

    # obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
    obj@se <- FindVariableFeatures(obj@se)
    obj@se <- RunPCA(obj@se, features = rownames(obj@se), seed.use = obj@params$rseed)
    obj@se <- RunUMAP(obj@se, features = rownames(obj@se), seed.use = obj@params$rseed)
    obj@se <- RunTSNE(obj@se, features = rownames(obj@se), seed.use = obj@params$rseed)
  }
  str(obj@se)


  return(obj)
})

# Exp1 ----
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
    if (FALSE) { # TESTING
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

    # if (!is.null(obj@params$pathw)) {
    #   # obj <- obj_updateParams(obj, updateCurrent = TRUE, pathw = pathw)
    #   obj <- obj_setGenes(obj)
    #   gene.flag <- rownames(obj@se) %in% obj@params$genes
    #   message("LOG: Number of genes: ", length(obj@params$genes))
    #   message("LOG: intersection pathw and genenames: ", sum(gene.flag))
    #   obj@se <- obj@se[gene.flag, ]

    #   obj@se <- SetAssayData(obj@se, layer = "scale.data", new.data = obj@se)
    #   # obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
    #   obj@se <- FindVariableFeatures(obj@se)
    #   obj@se <- RunPCA(obj@se, features = rownames(obj@se))
    #   obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
    # }

    message("Completed Loading")
    return(obj)
  }
)


# Exp2 -------
# Define the 'Exp2' Class that inherits from 'database'
setClass("Exp2",
  contains = "Exp"
)

## obj_loadData -------
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

    # if (!is.null(obj@params$pathw)) {
    #   obj <- obj_setGenes(obj)
    #   gene.flag <- rownames(obj@se) %in% obj@params$genes
    #   message("LOG: Number of genes: ", length(obj@params$genes))
    #   message("LOG: intersection pathw and genenames: ", sum(gene.flag))
    #   obj@se <- obj@se[gene.flag, ]

    #   obj@se <- SetAssayData(obj@se, layer = "scale.data", new.data = obj@se)
    #   # obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
    #   obj@se <- FindVariableFeatures(obj@se)
    #   obj@se <- RunPCA(obj@se, features = rownames(obj@se))
    #   obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
    # }

    message("Completed Loading")
    return(obj)
  }
)


# Exp3 ---------
setClass("Exp3",
  contains = "Exp"
)

## obj_loadData ----
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

    # if (!is.null(obj@params$pathw)) {
    #   obj <- obj_setGenes(obj)
    #   gene.flag <- rownames(obj@se) %in% obj@params$genes
    #   message("LOG: Number of genes: ", length(obj@params$genes))
    #   message("LOG: intersection pathw and genenames: ", sum(gene.flag))
    #   obj@se <- obj@se[gene.flag, ]

    #   obj@se <- SetAssayData(obj@se, layer = "scale.data", new.data = obj@se)
    #   # obj@se <- ScaleData(obj@se, features = rownames(obj@se), layer = "counts")
    #   obj@se <- FindVariableFeatures(obj@se)
    #   obj@se <- RunPCA(obj@se, features = rownames(obj@se))
    #   obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
    # }

    message("Completed Loading")
    return(obj)
  }
)



# obj_plotGoldUmap -----
setMethod("obj_plotGoldUmap", "Exp", function(obj) {
  umap_celltypes <- DimPlot(obj@se, reduction = "umap", group.by = "ctype")
  return(list("umap_celltypes" = umap_celltypes))
})
