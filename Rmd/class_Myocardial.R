# Myocardial ----
setClass("Myocardial",
  contains = "database"
)

## obj_loadData ----
setMethod(
  "obj_loadData", "Myocardial",
  function(obj,
             data_path = "/app/data/MyocardialInfarction/e61af320-303a-4029-8500-db6636bba0d4.rds",
           ...) {
    if (debug) message("DEBUG: obj_loadData |init load myocardial pathw is ", obj@params$pathw)
    if(FALSE){
      readRDS(data_path)
    }

    if (!is.null(obj@se.org)) {
      message("LOG: obj_loadData | Copy from original se")
      obj@se <- obj@se.org
    } else {
      message("LOG: obj_loadData | Full loading")

      obj <- obj_updateParams(obj, updateCurrent = TRUE, data_path = data_path)
      obj@se <- readRDS(data_path)
      str(obj@se)
      # message("\n\n\n\n")
      # # Show structure of each element in the new environment
      # for (element in names(new_env)) {
      #   cat(paste("Structure of element:", element, "\n"))
      #   str(new_env[[element]])
      #   cat("\n")
      # }
      
      # obj@se <- UpdateSeuratObject(MyocardialCortex)
      # rm(MyocardialCortex)
      # message("LOG: obj_loadData | updated Myocardial")
      options(future.globals.maxSize= 6*1024^3)
      obj@se$ctype <- Idents(obj@se)
      obj@se <- ScaleData(obj@se)
      obj@se <- FindVariableFeatures(obj@se)
      obj@se <- RunPCA(obj@se, features = rownames(obj@se))
      obj@se <- RunUMAP(obj@se, features = rownames(obj@se))
      str(obj@s)


      if (obj@params$test) {
        if (debug) message("DEBUG: obj_loadData | TEST selected -> reducing dataset")
        if (!is.null(obj@params$pathw)) {
          tgenes <- nrow(obj@se)
        } else {
          tgenes <- min(obj@params$test_genes, nrow(obj@se))
        }
        tsamples <- min(obj@params$test_samples, ncol(obj@se))

        metadata <- metadata[1:tsamples, ]

        obj@se <- obj@se[1:tgenes, 1:tsamples]
        obj@se <- obj@se[Matrix::rowSums(obj@se) > 0, Matrix::colSums(obj@se) > 0]
      }

      # save obj@se
      # obj@se.org <- obj@se
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
setMethod("obj_plotGoldUmap", "Myocardial", function(obj) {
  umap_celltypes <- DimPlot(obj@se, reduction = "umap", group.by = "ident")
  return(list(umap_celltypes = umap_celltypes))
})
