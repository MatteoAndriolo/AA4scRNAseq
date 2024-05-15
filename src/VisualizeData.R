source("src/loaddata.R")

createPlots <- function(filename, load.func){
  se=load.func()
  message("Dataset loaded")
  
  se = RunPCA(se, features = Seurat::VariableFeatures(se))
  message("RunPCA finished")

  #imgname=sprintf("%s/PCA.png",filename)
  # png(file = imgname, width = 600, height = 350)
  #message(sprintf("Saving Imagine --- %s", imgname))
  pdfname=sprintf("%s/PCA.pdf",filename)
  pdf(file = pdfname)
  message(sprintf("Saving PDF --- %s", pdfname))
  #PCAPlot(se)
  tryCatch({
    PCAPlot(se)
    #message(pca_result)
  }, error=function(e) {
      message("Failed to generate PCA plot: ", e$message)
  })
  dev.off()
  if (file.exists(imgname)) {
    print("File exists.")
  } else {
    print("File does not exist.")
  }

  se = RunUMAP(se, features = Seurat::VariableFeatures(se))
  message("RunUMAP finished")
  # imgname=sprintf("%s/UMAP.png",filename)
  # message(sprintf("Saving Imagine --- %s", imgname))
  # png(file = imgname, width = 600, height = 350)
  pdfname=sprintf("%s/PCA.pdf",filename)
  pdf(file = pdfname)
  message(sprintf("Saving PDF --- %s", pdfname))
  #pdf(file = imgname, width = 600, height = 350)
  tryCatch({
      UMAPPlot(se)
      #message(umap_result)
  }, error=function(e) {
      message("Failed to plot UMAP: ", e$message)
  })
  dev.off()
  if (file.exists(imgname)) {
    print("File exists.")
  } else {
    print("File does not exist.")
  }
  
  return(se)
}
  
#%% melanoma
filename="./data/Melanoma"
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


