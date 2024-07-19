# Mouse ----
setClass("Mouse",
  contains = "database"
)

## obj_loadData ----
setMethod(
  "obj_loadData", "Mouse",
  function(obj,
             data_path = "/app/data/MouseCortex/MouseCortex.RData",
           ...) {
    if (debug) message("DEBUG: obj_loadData |init load mouse pathw is ", obj@params$pathw)
    if(FALSE){
      load(data_path)
    }

    if (!is.null(obj@se.org)) {
      message("LOG: obj_loadData | Copy from original se")
      obj@se <- obj@se.org
    } else {
      message("LOG: obj_loadData | Full loading")

      obj <- obj_updateParams(obj, updateCurrent = TRUE, data_path = data_path)
      load(data_path)
      head(rownames(MouseCortex))
      message("LOG: obj_loadData | updating MouseCortex")
      
      MouseCortex <- UpdateSeuratObject(MouseCortex)
      head(rownames(MouseCortex))

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

      # save obj@se
      obj@se.org <- obj@se
    }

    if (!is.null(obj@params$pathw)) {
      if (debug) {
        message("DEBUG: obj_loadData | pathw is ", obj@params$pathw)
        message("DEBUG: obj_loadData | type of pathw ", typeof(obj@params$pathw))
      }
      message("LOG: obj_loadData | Loading pathw ", obj@params$pathw)
      obj <- obj_setGenes(obj)
      gene_names <- rownames(obj@se)
      message("gene names --")
      gene_names[1:5]
      message("gene names --")
      gene_names[1:5]
      gene.flag <- gene_names %in% obj@params$genes

      if (debug) {
        message("DEBUG: obj_loadData | pathw | number genes pathw = ", length(obj@params$genes))
        message("DEBUG: obj_loadData | pathw | first 5 ", toString(obj@params$genes[1:5]))
        message("DEBUG: obj_loadData | pathw | number genes data = ", length(gene_names))
        message("DEBUG: obj_loadData | pathw | first 5 ", toString(gene_names[1:5]))
        message("DEBUG: obj_loadData | pathw | intersectoin gives ", sum(gene.flag), " genes")
      }

      obj@se <- obj@se[gene.flag, ]
      message("LOG: obj_loadData | pathw | new dimension of se is ", dim(obj@se)[[1]], " ", dim(obj@se)[[2]])
    }

    message("LOG: obj_loadData | Completed Loading")
    return(obj)
  }
)

## obj_plotGoldUmap
setMethod("obj_plotGoldUmap", "Mouse", function(obj) {
  umap_celltypes <- DimPlot(obj@se, reduction = "umap", group.by = "idents")
  return(list(umap_celltypes = umap_celltypes))
})
