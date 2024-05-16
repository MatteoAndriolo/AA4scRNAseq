source("src/loaddata.R")
library(ggplot2)

createPlots <- function(filename, load.func){
  se=load.func()
  message("Dataset loaded")
  
  se = RunPCA(se, features = Seurat::VariableFeatures(se))
  message("RunPCA finished")
  
  message("Saving Imagine --- /app/elbow.pdf")
  #pdf(file = "/app/elbow.pdf")
  elbowplot=ElbowPlot(se)
  ggsave("/app/elbow.pdf", plot=elbowplot)
  #dev.off()



  imgname="/app/data/Melanoma/PCA.png"

  message(sprintf("Saving Imagine --- %s", imgname))
  
  sink(stdout(), type = "message")
  pcaplot=PCAPlot(se)
  ggsave(imgname, plot=pcaplot)
  #DimPlot(se, reduction = "pca")
  #dev.off()
  
  if (file.exists(imgname)) {
    print("File exists.")
  } else {
    print("File does not exist.")
    stop()
  }

  se = RunUMAP(se, features = Seurat::VariableFeatures(se))
  message("RunUMAP finished")
  imgname="/app/data/Melanoma/UMAP.png"
  message(sprintf("Saving Imagine --- %s", imgname))
  
  #png(file = imgname, width = 600, height = 350)
  umapplot = UMAPPlot(se)
  ggsave(imgname, plot=umapplot)
  #DimPlot(se, reduction = "umap")
  #dev.off()
  
  if (file.exists(imgname)) {
    print("File exists.")
  } else {
    print("File does not exist.")
  }
  
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

##%% Exp1
#filename="./data/AllonKleinLab/Experiment1"
#createPlots(filename,loadExp1)
#
##%% Exp2
#filename="./data/AllonKleinLab/Experiment2"
#createPlots(filename,loadExp2)
#
##%% Exp3
#filename="./data/AllonKleinLab/Experiment3"
#createPlots(filename,loadExp3)


