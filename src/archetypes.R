source("src/loaddata.R")

# ------------------------------------------------------------------------------
archetypesMelanoma <- function(k, test=FALSE){
  namefile = paste("out/Myocardite/Archetypes_",k,sep="")
  con <- file(paste(namefile,".Rlog",sep = ""))
  sink(con, append=TRUE)
  sink(con, appent)
  
              

  se = loadMelanoma()
  m=as.matrix(se@assays$scRNA@layers$counts)
  s=timestamp()
  a=archetypes::archetypes(m,k=k,maxIterations = 1)
  e=timestamp()
  # archetypes::xyplot(a,m)
  save(a,file=paste(namefile,".rds",sep=""))
  
}

archetypesMyocardialInfarction <- function(k, test=FALSE){
  namefile = paste("/app/out/MyocardialInfarction/Archetypes_",k,sep="")
  con=file(paste(namefile,".Rlog",sep=""))
  sink(con, append=TRUE)
  sink(con, append=TRUE, type="message")
  cat("Writing into sinkfile?",namefile,"\n")
  message("Writing message to output\n")
  #se=loadMyocardialInfarction()
  #m=as.matrix(se@assays$RNA@data )
  cat("MyocardialInfarction\n")
  s=timestamp()
  #a=archetypes::archetypes(m, k)
  e=timestamp()
  sink()
}

archetypesMouseCortex <- function(k, test=FALSE){
  namefile = paste("out/MouseCortex/Archetypes_",k,sep="")
  sinkToFile(paste(namefile,".log",sep=""))
  cat("Not Implemented")
  sinkDefault()
}

archetypesExp1 <- function(k, test=FALSE){
  namefile = paste("out/AllonKleinLab/Exp1/Archetypes_",k,sep="")
  sinkToFile(paste(namefile,".log",sep=""))
  cat("Not Implemented")
  sinkDefault()
}

archetypesExp2 <- function(k, test=FALSE){
  namefile = paste("out/AllonKleinLab/Exp2/Archetypes_",k,sep="")
  sinkToFile(paste(namefile,".log",sep=""))
  cat("Not Implemented")
  sinkDefault()
}

archetypesExp3 <- function(k, test=FALSE){
  namefile = paste("out/AllonKleinLab/Exp3/Archetypes_",k,sep="")
  sinkToFile(paste(namefile,".log",sep=""))
  cat("Not Implemented")
  sinkDefault()
}

# Function to run Archetype analysis
runArchetypeAnalysis <- function(dataset, numArchetypes) {
  print(paste("Running Archetype analysis on", dataset, "with", numArchetypes, "archetypes"))

  if(dataset=="Melanoma"){
    archetypesMelanoma(numArchetypes)
  }else if(dataset=="Myocardial"){
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

