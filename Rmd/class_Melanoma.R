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


      if (obj@params$test) {
        if (!is.null(obj@params$pathw)) {
          tgenes <- nrow(se)
        } else {
          tgenes <- min(obj@params$test_genes, nrow(se))
        }
        tsamples <- min(obj@params$test_samples, ncol(se))

        metadata <- metadata[1:tsamples, ]

        se <- se[1:tgenes, 1:tsamples]
        se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
      }

      obj@se <- CreateSeuratObject(counts = se, meta.data = metadata)
      if (debug) message("DEBUG: Seurat object has dimension ", dim(se)[[1]], " ", dim(se)[[2]])
      obj@se <- ScaleData(obj@se, layer = "counts")
      if (debug) message("DEBUG: Data scaled")
      obj@se <- FindVariableFeatures(obj@se)
      if (debug) message("DEBUG: Variable features found")
      obj@se <- RunPCA(obj@se, features = VariableFeatures(obj@se))
      if (debug) message("DEBUG: PCA done")
      obj@se <- RunUMAP(obj@se, features = VariableFeatures(obj@se))
      if (debug) message("DEBUG: UMAP done")
      obj@se.org <- obj@se
      if (debug) message("DEBUG: Seurat object created")

      if (obj@params$hvf) {
        message("LOG: HVF")
        obj@se <- obj@se[which(obj@se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
        obj@se <- obj@se[Matrix::rowSums(obj@se) > 0, Matrix::colSums(obj@se) > 0]
        message("LOG: HVF: new dimension of se is ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])
      }


      obj <- obj_updateParams(obj,
        updateCurrent = TRUE,
        data_path = data_path
      )
    } else {
      message("LOG: Copy from original se")
      obj@se <- obj@se.org
    }

    if (debug) message("DEBUG: pathw is ", obj@params$pathw)
    if (debug) message("DEBUG: type of pathw", typeof(obj@params$pathw))

    if (!is.null(obj@params$pathw)) {
      message("LOG: Loading pathw ", obj@params$pathw)
      obj <- obj_updateParams(obj, updateCurrent = TRUE, pathw = obj@params$pathw)
      obj <- obj_setGenes(obj)
      message("LOG: Number of genes: ", length(obj@params$genes))
      gene_names <- rownames(obj@se)
      gene.flag <- gene_names %in% obj@params$genes
      message("LOG: intersection pathw and seurat: ", sum(gene.flag))
      obj@se <- obj@se[obj@params$genes, ]
      message("LOG: dimension of se is ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])
    }

    message("LOG: Completed Loading")
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
  return(obj)
})
