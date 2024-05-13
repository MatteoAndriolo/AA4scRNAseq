source("src/loaddata.R")
setwd("/app")

createPlots <- function(filename, load.func){
  message(sprintf("Loading %s ", filename))
  se=load.func()
  message("Dataset loaded")
  
  se = RunPCA(se, features = Seurat::VariableFeatures(se))
  message("RunPCA finished")
  se = RunUMAP(se, features = Seurat::VariableFeatures(se))
  message("RunUMAP finished")
  
  message(sprintf("writing to %s/PCA.png", filename))
  png(file = sprintf("%s/PCA.png",filename), width = 600, height = 350)
  PCAPlot(se)
  dev.off()
  
  png(file = sprintf("%s/UMAP.png",filename), width = 600, height = 350)
  UMAPPlot(se)
  dev.off()

  message(sprintf("Finished %s", filename))
}
  
# #%% melanoma
# filename="/app/data/Melanoma"
# createPlots(filename, loadMelanoma)

# #%% Miocardial
# filename="/app/data/MyocardalInarction"
# createPlots(filename, loadMyocardialInfarction)

# #%% MouseCortex
# filename="/app/data/MouseCortex"
# createPlots(filename, loadMyocardialInfarction)

#%% Exp1
# filename="/app/data/AllonKleinLab/Experiment1"
# createPlots(filename,loadExp1)
# 
# #%% Exp2
# filename="/app/data/AllonKleinLab/Experiment2"
# createPlots(filename,loadExp2)

#%% Exp3
filename="/home/andrioloma/MasterThesis/data/AllonKleinLab/Experiment3"
createPlots(filename,loadExp3)


