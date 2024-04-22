source("src/loaddata.R")

archetypesMelanoma <- function(){
  se = loadMelanoma()
  m=as.matrix(se@assays$scRNA@layers$counts)
  s=timestamp()
  a=archetypes::archetypes(m,k=5,maxIterations = 1)
  e=timestamp()
  archetypes::xyplot(a,m)
  #return(se)
}


se=archetypesMelanoma()
