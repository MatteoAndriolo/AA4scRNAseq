# Melanoma ----
setClass("Melanoma",
  contains = "database"
)

## obj_loadData ----
setMethod(
  "obj_loadData", "Melanoma",
  function(obj,
           data_path = "/app/data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt",
           ...) {
    if (debug) message("DEBUG: obj_loadData |INITLOADMELANOMA pathw is ", obj@params$pathw)

    if (!is.null(obj@se.org)) {
      message("LOG: obj_loadData | Copy from original se")
      obj@se <- obj@se.org
    } else {
      message("LOG: obj_loadData | Full loading")
      obj <- obj_updateParams(obj, updateCurrent = TRUE, data_path = data_path)
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
      if (debug) message("DEBUG: obj_loadData | data matrix has dimension ", dim(se)[[1]], " ", dim(se)[[2]])
      se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
      if (debug) message("DEBUG: obj_loadData | data matrix has dimension post rem0 ", dim(se)[[1]], " ", dim(se)[[2]])


      if (obj@params$test) {
        if (debug) message("DEBUG: obj_loadData | TEST selected -> reducing dataset")
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
      if (debug) message("DEBUG: obj_loadData | Seurat object has dimension ", dim(se)[[1]], " ", dim(se)[[2]])
      obj@se <- ScaleData(obj@se, layer = "counts")
      obj@se <- FindVariableFeatures(obj@se)
      obj@se <- RunPCA(obj@se, features = rownames(obj@se))
      obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
      # obj@se <- RunPCA(obj@se, features = VariableFeatures(obj@se))
      # obj@se <- RunUMAP(obj@se, features = VariableFeatures(obj@se))

      # save obj@se
      obj@se.org <- obj@se

      # HVF
      if (obj@params$hvf) {
        message("LOG: obj_loadData | HVF")
        obj@se <- obj@se[which(obj@se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
        obj@se <- obj@se[Matrix::rowSums(obj@se) > 0, Matrix::colSums(obj@se) > 0]
        message("LOG: obj_loadData | HVF: new dimension of se is ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])

        # message("LOG: obj_loadData | obj_loadData | rescale, hvf, reduce after hvf")
        # obj@se <- ScaleData(obj@se, layer = "counts")
        # obj@se <- FindVariableFeatures(obj@se, features = rownames(obj@se))
        obj@se <- RunPCA(obj@se, features = rownames(obj@se))
        obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
      }
    }

    if (!is.null(obj@params$pathw)) {
      if (debug) {
        message("DEBUG: obj_loadData | pathw is ", obj@params$pathw)
        message("DEBUG: obj_loadData | type of pathw ", typeof(obj@params$pathw))
      }
      message("LOG: obj_loadData | Loading pathw ", obj@params$pathw)
      obj <- obj_setGenes(obj)
      gene_names <- rownames(obj@se)
      gene.flag <- gene_names %in% obj@params$genes

      if(debug){
        message("DEBUG: obj_loadData | pathw | number genes pathw = ", length(obj@params$genes) )
        message("DEBUG: obj_loadData | pathw | first 5 ", toString(obj@params$genes[1:5]))
        message("DEBUG: obj_loadData | pathw | number genes data = ", length(gene_names ))
        message("DEBUG: obj_loadData | pathw | first 5 ", toString(gene_names[1:5]))
        message("DEBUG: obj_loadData | pathw | intersectoin gives ", sum(gene.flag), " genes")
      }

      obj@se <- obj@se[gene.flag , ]
      message("LOG: obj_loadData | pathw | new dimension of se is ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])

      # message("LOG: obj_loadData | rescale, hvf, reduce after pathw")
      # obj@se <- ScaleData(obj@se, layer = "counts")
      # obj@se <- FindVariableFeatures(obj@se)
      # obj@se <- RunPCA(obj@se, features = rownames(obj@se))
      # obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
      # obj@se <- RunPCA(obj@se, features = VariableFeatures(obj@se))
      # obj@se <- RunUMAP(obj@se, features = VariableFeatures(obj@se))
    }

    message("LOG: obj_loadData | Completed Loading")
    return(obj)
  }
)

## obj_getSeData ----
setMethod("obj_getSeData", "Melanoma", function(obj) {
  if (debug) message("DEBUG: obj_getSeData | return counts")
  return(obj@se@assays$RNA@layers$counts)
})

## obj_plotGoldUmap
setMethod("obj_plotGoldUmap", "Melanoma", function(obj) {
  ct <- "non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK."
  umap_celltypes <- DimPlot(obj@se, reduction = "umap", group.by = ct)
  umap_tumor <- DimPlot(obj@se, reduction = "umap", group.by = "tumor")
  return(list(umap_celltypes = umap_celltypes, umap_tumor = umap_tumor))
})

### obj_getCellTypesList ----
# setMethod("obj_getCellTypesList", "Melanoma", function(obj) {
#  ct <- "non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK."
#  return(obj@se[[ct]])
# })
#
### obj_getCellTypesMetaDataName ----
# setMethod("obj_getCellTypesMetaDataName", "Melanoma", function(obj) {
#  ct <- "non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK."
#  return(ct)
# })
#
