# rstudioapi::filesPaneNavigate("/app")
# setwd("/app")
source("src/loaddata.R")
# ------------------------------------------------------------------------------
archetypesMelanoma <- function(k, small=TRUE){
  se = loadMelanoma()
  
  if(small){
    m = se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank>0),]
  }else{
    m=se@assays$scRNA@layers$counts
  }
  m=as.matrix(m)
  
  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
  e=timestamp()
  
  
  save(a,file=paste("/app/out/Melanoma/Archetypes_",k,".rds",sep=""))
}

# ------------------------------------------------------------------------------
archetypesMyocardialInfarction <- function(k, small=TRUE){
  se=loadMyocardialInfarction()
  m=as.matrix(se@assays$RNA@data)

  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
  e=timestamp()

  save(a,file=paste("/app/out/MyocardialInfarction/Archetypes_",k,".rds",sep=""))
}

# ------------------------------------------------------------------------------
archetypesMouseCortex <- function(k, small=TRUE){
  se=loadMouseCortex()
  # se is Seurat4 object
  warning("Error problem with negative data")
  stop()
#  m=as.matrix(se@data)
#  
#  s=timestamp()
#  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
#  e=timestamp()
#  
#  cat("Not Implemented")
#  sinkDefault()
}

# ------------------------------------------------------------------------------
archetypesExp1 <- function(k,small=TRUE){
  se=loadExp1()
  
  if(small){
    m = se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank>0),]
  }else{
    m = se@assays$RNA@layers$counts
  }
  m=as.matrix(m)
  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
  e=timestamp()
  
  save(a,file=paste("/app/out/AllonKleinLab/Experiment1/Archetypes_",k,".rds",sep=""))
}

# ------------------------------------------------------------------------------
archetypesExp2 <- function(k, small=TRUE){
  se=loadExp2()
  
  if(small){
    m = se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank>0),]
  }else{
    m = se@assays$RNA@layers$counts
  }
  m=as.matrix(m)
  
  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
  e=timestamp()
  save(a,file=paste("/app/out/AllonKleinLab/Experiment2/Archetypes_",k,".rds",sep=""))
}

# ------------------------------------------------------------------------------
archetypesExp3 <- function(k,small=TRUE){
  se=loadExp3()
  if(small){
    m = se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank>0),]
  }else{
    m = se@assays$RNA@layers$counts
  }
  m=as.matrix(m)
  
  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=1)
  e=timestamp()
  save(a,file=paste("/app/out/AllonKleinLab/Experiment3/Archetypes_",k,".rds",sep=""))
}

# ------------------------------------------------------------------------------
# Function to run Archetype analysis
# ------------------------------------------------------------------------------
runArchetypeAnalysis <- function(dataset, numArchetypes) {
  print(paste("Running Archetype analysis on", dataset, "with", numArchetypes, "archetypes"))

  if(dataset=="Melanoma"){
    archetypesMelanoma(numArchetypes)
  }else if(dataset=="MyocardialInfarction"){
    archetypesMyocardialInfarction(numArchetypes)
  }else if(dataset=="MouseCortex"){
    archetypesMouseCortex(numArchetypes)
  }else if(dataset=="AllonKleinLab/Experiment1"){
    archetypesExp1(numArchetypes)
  }else if(dataset=="AllonKleinLab/Experiment2"){
    archetypesExp2(numArchetypes)
  }else if(dataset=="AllonKleinLab/Experiment3"){
    archetypesExp3(numArchetypes)
  }else{
    cat("Dataset is Wrong")
  }
}

