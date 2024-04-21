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

loadExp1 <- function(){
  requite(Matrix)
  require(dplyr)
  clone_matrix <- Matrix::readMM("data/AllonKleinLab/Experiment1/stateFate_inVitro_clone_matrix.mtx")
  normed_count <- Matrix::readMM("data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx")
  metadata <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt",header=TRUE,sep = "\t" )
  gene_names <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")
  neutrophil_monocyte_trajectory <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_monocyte_trajectory.txt",header=TRUE)
  neutrophil_pseudotime <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_pseudotime.txt",header=TRUE)
}

se=loadMelanoma()
