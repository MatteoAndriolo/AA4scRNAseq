if(getwd()!="/app"){
  setwd("/app")
}
getwd()
source("src/loaddata.R")

con <- file("test.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

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


