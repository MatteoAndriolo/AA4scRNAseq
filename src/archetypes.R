source("src/loaddata.R")

# ------------------------------------------------------------------------------
archetypesMelanoma <- function(k){
  se = loadMelanoma()
  m=as.matrix(se@assays$scRNA@layers$counts)

  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
  e=timestamp()

  save(a,file=paste("/app/out/Melanoma/Archetypes_",k,".rds",sep=""))
}

# ------------------------------------------------------------------------------
archetypesMyocardialInfarction <- function(k){
  se=loadMyocardialInfarction()
  m=as.matrix(se@assays$RNA@data)
  rm(se)
  gc()

  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
  e=timestamp()

  save(a,file=paste("/app/out/MyocardialInfarction/Archetypes_",k,".rds",sep=""))
}

# ------------------------------------------------------------------------------
archetypesMouseCortex <- function(k){
  se=loadMouseCortex()
  # se is Seurat4 object
  m=as.matrix(se@data)
  
  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
  e=timestamp()
  
  cat("Not Implemented")
#  sinkDefault()
}

# ------------------------------------------------------------------------------
archetypesExp1 <- function(k){
  m=loadExp1()
  m=as.matrix(m)
  
  s=timestamp()
  a=archetypes::archetypes(m, k=k, verbose=TRUE, maxIterations=10)
  e=timestamp()
  
  save(a,file=paste("/app/out/AllonKleinLab/Experiment1/Archetypes_",k,".rds",sep=""))
}

# ------------------------------------------------------------------------------
archetypesExp2 <- function(k){
  se=loadExp2()
}

# ------------------------------------------------------------------------------
archetypesExp3 <- function(k){
  se=loadExp3()
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
  }else if(dataset=="Exp1"){
    archetypesExp1(numArchetypes)
  }else if(dataset=="Exp2"){
    archetypesExp2(numArchetypes)
  }else if(dataset=="Exp3"){
    archetypesExp3(numArchetypes)
  }else{
    cat("Dataset is Wrong")
  }
}

