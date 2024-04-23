# SMALL SNAPSHOT
loadMelanoma <- function(){
  library(dplyr)
  library(Matrix)
  
  a <- read.table("data/MelanomaSomethingPiccolo/GSE72056_melanoma_single_cell_revised_v2.txt", header = TRUE)
  a=a[!duplicated(a[,1]),]
  
  rownames(a)=a[,1]
  a=a[,2:ncol(a)]
  
  metadata=a[1:3,]
  metadata <- t(metadata)%>% 
    data.frame() %>% 
    mutate(across(where(is.character), as.numeric)) 
  
  ge=a[4:nrow(a),]
  ge <- ge %>% 
    data.frame() %>% 
    mutate(across(where(is.character), as.numeric)) 
  rm(a)
  
  se=SeuratObject::CreateSeuratObject(ge,assay="scRNA",project="MelanomaSC",meta.data=metadata)
  return(se)
}


# LARGE SNAPSHOT
loadMyocardialInfarction <- function(){
  library(Matrix)
  se=readRDS("data/MyocardialInfarction/e61af320-303a-4029-8500-db6636bba0d4.rds")
  return(se)
}

# SMALL TIMESERIES
loadMouseCortex <- function(){
  load("data/MouseCortex/MouseCortex.RData")
  MouseCortex=MouseCortex
  return(MouseCortex)
}

# LARGE TIMESERIES
loadExp1 <- function(){
  requite(Matrix)
  require(dplyr)
  require(monocle3)
  clone_matrix <- Matrix::readMM("data/AllonKleinLab/Experiment1/stateFate_inVitro_clone_matrix.mtx")
  normed_count <- Matrix::readMM("data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx")
  cell_metadata <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt",header=TRUE,sep = "\t" )
  gene_names <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")
  neutrophil_monocyte_trajectory <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_monocyte_trajectory.txt",header=TRUE,sep="\t")
  neutrophil_pseudotime <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_pseudotime.txt",header=TRUE, sep="\t")
  
  # TODO create Seurat / Monocle3 + SueratWrapper
  # as in https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/
}

loadExp2 <- function(){
  clone_matrix <- Matrix::readMM("data/AllonKleinLab/Experiment2/stateFate_inVivo_clone_matrix.mtx")
  gene_names <- read.table("data/AllonKleinLab/Experiment2/stateFate_inVivo_gene_names.txt",sep="\t")
  metadata <- read.table("data/AllonKleinLab/Experiment2/stateFate_inVivo_metadata.txt",sep="\t")
  normed_counts <- Matrix::readMM("data/AllonKleinLab/Experiment2/stateFate_inVivo_normed_counts.mtx")
  
}
loadExp3 <- function(){
  clone_matrix = Matrix::readMM("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_clone_matrix.mtx")
  gene_names <- read.table( "data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_gene_names.txt")
  metadata <- read.table("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_metadata.txt", sep="\t")
  normed_counts <- Matrix::readMM("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_normed_counts.mtx")
}
