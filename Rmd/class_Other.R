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
