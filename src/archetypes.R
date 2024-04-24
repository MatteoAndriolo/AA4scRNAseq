source("src/loaddata.R")

IF_ARCHETYPES <- function(matrix, k, maxiter=1){
  return(
}
# ------------------------------------------------------------------------------
archetypesMelanoma <- function(k){
  se = loadMelanoma()
  m=as.matrix(se@assays$scRNA@layers$counts)

  s=timestamp()
  a=archetypes::archetypes(matrix, k=k, verbose=TRUE, maxIterations=1)
  e=timestamp()

  save(a,file=paste("/app/out/Melanoma/Archetypes_",k,".rds",sep=""))
}

archetypesMyocardialInfarction <- function(k){
  se=loadMyocardialInfarction()
  m=as.matrix(se@assays$RNA@data)
  rm(se)
  gc()

  s=timestamp()
  a=archetypes::archetypes(matrix, k=k, verbose=TRUE, maxIterations=1)
  e=timestamp()

  save(a,file=paste("/app/out/MyocardialInfarction/Archetypes_",k,".rds",sep=""))
}

archetypesMouseCortex <- function(k){
#  namefile = paste("out/MouseCortex/Archetypes_",k,sep="")
#  sinkToFile(paste(namefile,".log",sep=""))
  cat("Not Implemented")
#  sinkDefault()
}

archetypesExp1 <- function(k){
#  namefile = paste("out/AllonKleinLab/Exp1/Archetypes_",k,sep="")
#  sinkToFile(paste(namefile,".log",sep=""))
  cat("Not Implemented")
#  sinkDefault()
}

archetypesExp2 <- function(k){
#  namefile = paste("out/AllonKleinLab/Exp2/Archetypes_",k,sep="")
#  sinkToFile(paste(namefile,".log",sep=""))
  cat("Not Implemented")
#  sinkDefault()
}

archetypesExp3 <- function(k){
#  namefile = paste("out/AllonKleinLab/Exp3/Archetypes_",k,sep="")
#  sinkToFile(paste(namefile,".log",sep=""))
  cat("Not Implemented")
#  sinkDefault()
}

# Function to run Archetype analysis
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

