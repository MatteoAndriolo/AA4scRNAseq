# Name: unique
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")
debug <- TRUE
params <- list()
params$classname <- Sys.getenv("CLASSNAME")
params$debug <- debug 
params$hvf <- as.logical(Sys.getenv("HVF", "FALSE"))
params$k <- as.numeric(Sys.getenv("K", "8"))
params$max_iterations <- as.numeric(Sys.getenv("MAX_ITERATIONS", "100"))
params$num_restarts <- as.numeric(Sys.getenv("NUM_RESTARTS", "10"))
params$nworkers <- as.numeric(Sys.getenv("NWORKERS", "10"))
params$out_path <- Sys.getenv("OUT_PATH")
params$pathw <- as.numeric(Sys.getenv("PATHW", unset = "-1"))
params$test <- as.logical(Sys.getenv("TEST", "FALSE"))
params$test_genes <- as.numeric(Sys.getenv("TEST_GENES", "300"))
params$test_samples <- as.numeric(Sys.getenv("TEST_SAMPLES", "500"))

# params$nworkers <- parallel::detectCores() - 2
plan("multicore", workers = params$nworkers)

################### FIXING PARAMETERS FOR PRESENTATION 
#TODO remove this
# FOR NOW WE WILL NOT USE PATHW BUT ONLY HVF!!!
params$pathw <- -1
if (params$pathw > 0) {
  params$pathw <- pathways[[params$pathw]]
} else {
  params$pathw <- NULL
}

#if(is.null(params$k) & classname=="Melanoma"){
params$kappas<- 7:9
#}

params$hvf <- TRUE
################### END FIXING PARAMETERS FOR PRESENTATION 

for (k in names(params)) {
  message("LOG: main | param ", k, ": ", params[[k]])
}

# CREATE OBJECT -----
obj <- new(params$classname)
obj <- do.call(obj_updateParams, c(list(obj = obj), params))

# LOADING DATA -----
message("LOG: main | Loading Data")
obj <- obj_loadData(obj)
message("LOG: main | Loading Data Done")

# if(!is.null(params$pathw)){
#   message("LOG: main | starting with PATHW ", params$pathw)
#   obj <- obj_updateParams(obj, pathw = params$pathw)
#   message("LOG: main | reloading data with pathw")
#   obj <- obj_loadData(obj)
# }

# PERFORM ARCHETYPES ---------------------------------
message("LOG: main | Performing Archetypes")
obj@params$kappas=7:9
obj <- obj_performArchetypes(obj, doparallel = FALSE)
message("LOG: main | Performing Archetypes Done")

# message("LOG: main | assign AA clusters")
# obj <- obj_assignAAClusters(obj)
# message("LOG: main | assign AA clusters Done")

# SEURAT CLUSTERIZATOIN -------------------------------
message("LOG: main | Seurat Clustering")
obj <- obj_seuratCluster(obj)
message("LOG: main | Seurat Clustering Done")

# # VISUALIZE -------------------------------------------
# # Visualize Dataset
# message("LOG: main | Visualizing Data")
# obj <- obj_visualizeData(obj)
# message("LOG: main | Visualizing Data Done")
# 
# # Visualize Archetypes
# message("LOG: main | Visualizing Archetypes")
# obj <- obj_visualizeArchetypes(obj)
# message("LOG: main | Visualizing Archetypes Done")
# 
# # Umap Archetypes Plot
# message("LOG: main | Umap Archetypes")
# obj <- obj_umapArchetypes(obj)
# message("LOG: main | Umap Archetypes Done")
# 
# 
# 
# obj@se@meta.data$seurat_clusters
# obj@plots$umap_orig.ident <- DimPlot(obj@se, reduction = "umap", group.by = "orig.ident")
# obj@plots$umap_tumor <- DimPlot(obj@se, reduction = "umap", group.by = "tumor")
# obj@plots$umap_seucl <- DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
# obj@plots$umap_aacl <- DimPlot(obj@se, reduction = "umap", group.by = "aa_clusters")
# 
# # Compare aa_clusters and seurat_clasters with Ident()
# obj@compare$aa.orig <- table(obj@se@meta.data$aa_clusters, obj@se@meta.data$orig.ident)
# obj@compare$se.orig <- table(obj@se@meta.data$seurat_clusters, obj@se@meta.data$orig.ident)
# obj@compare$aa.se <- table(obj@se@meta.data$aa_clusters, obj@se@meta.data$seurat_clusters)

message("LOG: main | Saving Object")
obj_saveObj(obj, name="unique")
message("LOG: main | Saving Object Done")
