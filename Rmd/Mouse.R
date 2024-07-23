library(Seurat)
data_path <- "/app/data/MouseCortex/MouseCortex.RData"
load(data_path)
se <- CreateSeuratObject(MouseCortex@raw.data, meta.data = MouseCortex@meta.data, min.cells = 0, min.features = 0, names.field = 1, names.delim = "_")
Idents(se) <- MouseCortex@ident
rm(MouseCortex)
se <- NormalizeData(se, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
se <- ScaleData(se, layer = "counts", model.use = "linear", use.umi = FALSE, do.scale = TRUE, do.center = TRUE, scale.max = 10, block.size = 1000, min.cells.to.block = 3000, assay = "RNA")
se <- FindVariableFeatures(se)
se <- RunPCA(se, features = rownames(se))
se <- RunUMAP(se, features = rownames(se))

# plot= DimPlot(se, reduction="umap")
# ggsave(plot=plot, filename="/app/RmdTest/umapplot.png")
