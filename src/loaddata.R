require(Matrix)
require(dplyr)
require(Seurat)
# require(monocle3)

# SMALL SNAPSHOT
loadMelanoma <- function(hvf.nfeatures = 2000){
  # -------- BASIC SETUP
  ge <- read.table("data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt", header = TRUE)
  ge=ge[!duplicated(ge[,1]),] # remove duplicated genes
  rownames(ge)=ge[,1]         # extract from matrix rownames
  ge=ge[,2:ncol(ge)]          # elide rownames from gene expression matrix
  
  metadata=ge[1:3,]          # extract metadata
  metadata <- t(metadata) %>%
    data.frame() %>% 
    mutate(across(where(is.character), as.numeric)) 
  
  ge=ge[4:nrow(ge),]
  ge <- ge %>% 
    data.frame() %>% 
    mutate(across(where(is.character), as.numeric)) 
  # -------- END BASIC SETUP
  
  # -------- SEURAT
  se=SeuratObject::CreateSeuratObject(ge,project="MelanomaSC",meta.data=metadata)
  # Skip normalization since data is already normalized
  se=Seurat::ScaleData(se, layer= "counts") #,check.for.norm=FALSE)
  se=Seurat::FindVariableFeatures(se, hvf.nfeatures = hvf.nfeatures)
  
  #se=Seurat::RunPCA(se, features = Seurat::VariableFeatures(se))
  #se = RunUMAP(se, features = Seurat::VariableFeatures(se))
  #Seurat::UMAPPlot(se)
  
  # se=Seurat::RunUMAP(se, dims=9:15)
  # Seurat::DimPlot(se, reduction="pca")
  # Seurat::DimHeatmap(se, dims=1, cells=1000, balanced=TRUE)
  # Seurat::DimHeatmap(se, dims = 1:15, cells = 1000, balanced = TRUE)
  
  #ElbowPlot(se)
  
  return(se)
}

# LARGE SNAPSHOT
loadMyocardialInfarction <- function(hvf.nfeatures = 2000){
  se=readRDS("data/MyocardialInfarction/e61af320-303a-4029-8500-db6636bba0d4.rds")
  se=FindVariableFeatures(se)
  warning("Implementation Not completed due to problem with data")
  stop()
  return(se)
}

# SMALL TIMESERIES
loadMouseCortex <- function(hvf.nfeatures=2000){
  load("data/MouseCortex/MouseCortex.RData")
  MouseCortex=MouseCortex
  warning("Implementation Not completed due to problem with data")
  stop()
}

# LARGE TIMESERIES
loadExp1 <- function(hvf.nfeatures = 2000){
  # Binary matrix indicating clonal membership of each cell
  # # The rows of this file represent cells and correspond to the rows of _counts_matrix_in_vitro_ (above).
  # # The columns represent clones. Not every cell belongs to a clone. 
  #clone_matrix <- Matrix::readMM("data/AllonKleinLab/Experiment1/stateFate_inVitro_clone_matrix.mtx")
  
  # Number of transcripts for each gene in each cell
  se <- Matrix::readMM("/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx" )
  se <- CreateSeuratObject(counts=se)
  se=ScaleData(se, layer= "counts") #,check.for.norm=FALSE)
  se=FindVariableFeatures(se)
  
  # # cell metadata : cell type annotation
  # cell_metadata <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt",header=TRUE,sep = "\t" )
  # 
  # gene_names <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")
  # 
  # # List of cells belonging to the neutrophil/monocyte trajectory that were used in becnmark analysis
  # neutrophil_monocyte_trajectory <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_monocyte_trajectory.txt",header=TRUE,sep="\t")
  
  # # pseudotime for neutrophil trajectory cells
  # neutrophil_pseudotime <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_pseudotime.txt",header=TRUE, sep="\t")
  
  return(se)
}

loadExp2 <- function(hvf.nfeatures=2000){
  # clone_matrix <- Matrix::readMM("data/AllonKleinLab/Experiment2/stateFate_inVivo_clone_matrix.mtx")
  # gene_names <- read.table("data/AllonKleinLab/Experiment2/stateFate_inVivo_gene_names.txt",sep="\t")
  # metadata <- read.table("data/AllonKleinLab/Experiment2/stateFate_inVivo_metadata.txt",sep="\t")
  
  se <- Matrix::readMM("data/AllonKleinLab/Experiment2/stateFate_inVivo_normed_counts.mtx")
  se <- as(se,"CsparseMatrix")
  se <- CreateSeuratObject(counts=se)
  se=FindVariableFeatures(se)
  
  return(se) 
}

loadExp3 <- function(hvf.nfeatures = 2000){
  # clone_matrix = Matrix::readMM("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_clone_matrix.mtx")
  # gene_names <- read.table( "data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_gene_names.txt")
  # metadata <- read.table("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_metadata.txt", sep="\t")
  
  se <- Matrix::readMM("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_normed_counts.mtx")
  se <- CreateSeuratObject(counts=se)
  se=Seurat::ScaleData(se, layer= "counts") #,check.for.norm=FALSE)
  se=FindVariableFeatures(se)
  
  return(se)
}
