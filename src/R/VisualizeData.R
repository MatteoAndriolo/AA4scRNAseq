source("src/loaddata.R")
library(ggplot2)

createPlots <- function(filename, load.func){
  se=load.func()
  message("Dataset loaded")
  

  message("RunPCA finished")

  imgname=sprintf("%s/PCA.png", filename)
  message(sprintf("Saving Imagine --- %s",imgname))
  pcaplot=PCAPlot(se)
  ggsave(imgname, plot=pcaplot)

  imgname=sprintf("%s/UMAP.png", filename)
  message(sprintf("Saving Imagine --- %s",imgname))
  umapplot = UMAPPlot(se)
  ggsave(imgname, plot=umapplot)
  
  imgname=sprintf("%s/elbow.pdf", filename)
  message(sprintf("Saving Imagine --- %s",imgname))
  elbowplot=ElbowPlot(se)
  ggsave(imgname, plot=elbowplot)
  
  return(se)
}
  
#%% melanoma
filename="/app/data/Melanoma"
createPlots(filename, loadMelanoma)

# #%% Miocardial
# filename="/app/data/MyocardalInarction"
# createPlots(filename, loadMyocardialInfarction)

# #%% MouseCortex
# filename="/app/data/MouseCortex"
# createPlots(filename, loadMyocardialInfarction)

#%% Exp1
filename="./data/AllonKleinLab/Experiment1"
createPlots(filename,loadExp1)

#%% Exp2
filename="./data/AllonKleinLab/Experiment2"
createPlots(filename,loadExp2)

#%% Exp3
filename="./data/AllonKleinLab/Experiment3"
createPlots(filename,loadExp3)


