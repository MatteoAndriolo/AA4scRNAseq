source("src/loaddata.R")
library(Seurat)
# ------------------------------------------------------------
# Standard Seurat Pipeline
# ------------------------------------------------------------
stdSeurat <- function(raw.data) {
  source("src/loaddata.R")
  library(Seurat)
  se <- loadMouseCortex()
  VlnPlot(se, features = c("nCount_scRNA", "nFeature_scRNA"))
  plot1 <- FeatureScatter(se, feature1 = "nCount_scRNA", feature2 = "nFeature_scRNA")
  plot1
  se <- subset(se, subset = nFeature_scRNA < 10000 && nFeature_scRNA > 200 && nCount_scRNA < 20000 && nCount_scRNA > 4000)
  VlnPlot(se, features = c("nCount_scRNA", "nFeature_scRNA"))
  se <- SCTransform(se)

  # gene expression
  # Cell clusters were identified based on the Louvain clustering.
  # Differentially expressed genes (DEGs) in each cluster were identified using the “FindAllMarkers” function with default parameters.
  # Cell-type annotation was performed similarly to previous reports.14, 15
  # from https://onlinelibrary.wiley.com/doi/10.1111/cpr.13201
}

# Function to run Seurat analysis
seuratMelanoma <- function() {
  se <- loadMelanoma(hvf.nfeatures = 2000)
  Seurat::ElbowPlot(se)

  se <- Seurat::FindNeighbors(se, dims = 6:15)
  se <- Seurat::FindClusters(se, resolution = 0.5)
  # head(Seurat::Idents(se), 5)

  se <- RunUMAP(se, dims = 9:15)
  DimPlot(se, reduction = "umap")
  # FeaturePlot(se,features=c("tumor","malignant"))
  # VlnPlot(se, features = "tumor")

  se[["RNA"]]$counts %>%
    as.data.frame() %>%
    filter(rowname %in% Seurat::VariableFeatures(se))
}

seuratMyocardialInfarction <- function() {
  se <- loadMyocardialInfarction()
  se <- FindVariableFeatures(se, nfeatures = 2000)
}
seuratMouseCortex <- function() {
  se <- loadMouseCortex()
}
seuratExp1 <- function() {

}
seuratExp2 <- function() {

}
seuratExp3 <- function() {

}



runSeuratAnalysis <- function(dataset) {
  print(paste("Running Seurat analysis on", dataset))

  if (dataset == "Melanoma") {
    seuratMelanoma()
  } else if (dataset == "MyocardialInfarction") {
    seuratMyocardialInfarction()
  } else if (dataset == "MouseCortex") {
    seuratMouseCortex()
  } else if (dataset == "Exp1") {
    seuratExp1()
  } else if (dataset == "Exp2") {
    seuratExp2()
  } else if (dataset == "Exp3") {
    seuratExp3()
  } else {
    cat("Dataset is Wrong")
  }
}
