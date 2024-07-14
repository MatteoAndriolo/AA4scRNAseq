# Name: unique
sink(stdout(), type = "message")
sink(stdout(), type = "output")
# Import necessary R scripts
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")
# Fetch params from environment variables
params <- list()
params$classname <- Sys.getenv("CLASSNAME")
params$debug <- as.logical(Sys.getenv("DEBUG", "FALSE"))
params$hvf <- as.logical(Sys.getenv("HVF", "FALSE"))
params$max_iterations <- as.numeric(Sys.getenv("MAX_ITERATIONS", "100"))
params$num_restarts <- as.numeric(Sys.getenv("NUM_RESTARTS", "10"))
params$method <- Sys.getenv("METHOD")
params$nworkers <- as.numeric(Sys.getenv("NWORKERS", "1"))
params$out_path <- Sys.getenv("OUT_PATH")
params$pathw <- as.numeric(Sys.getenv("PATHW", "0"))
params$test <- as.logical(Sys.getenv("TEST", "FALSE"))
params$test_genes <- as.numeric(Sys.getenv("TEST_GENES", "300"))
params$test_samples <- as.numeric(Sys.getenv("TEST_SAMPLES", "500"))

# params$name <- "unique"
# params$k <- as.numeric(Sys.getenv("K", "8"))
# params$nworkers <- parallel::detectCores() - 2
debug <- params$debug
params$path_figures <- file.path(params$out_path, "figures")
plan("multicore", workers = params$nworkers)
if (params$pathw == 0) {
  params$pathw <- NULL
}

################### FIXING PARAMETERS FOR PRESENTATION
# TODO remove this
# FOR NOW WE WILL NOT USE PATHW BUT ONLY HVF!!!
# params$pathw <- 0
# if (params$pathw > 0) {
#  params$pathw <- pathways[[params$pathw]]
# } else {
#  params$pathw <- NULL
# }
#
# if(is.null(params$k) & classname=="Melanoma"){
params$kappas <- 4:8
# }
#
# params$hvf <- TRUE
# params$name <- "5.13.FS.unique."
# params$doFurthestSum <- TRUE
################### END FIXING PARAMETERS FOR PRESENTATION

# CREATE OBJECT -----
obj <- new(params$classname)
obj <- do.call(obj_updateParams, c(list(obj = obj), params))

# LOADING DATA -----
message("LOG: main | Loading Data")
obj <- obj_loadData(obj)
if (debug) message("DEBUG: main | Loading Data Done")

# PERFORM ARCHETYPES ---------------------------------
## ARCHETYPES 
if (obj@params$method == "archetypes"){
message("LOG: main | Performing Archetypes")
obj <- obj_performArchetypes(obj, doparallel = FALSE)
message("LOG: main | Performing Archetypes Done")

message("LOG: main | assign AA clusters")
obj <- obj_assignArchetypesClusters(obj)
message("LOG: main | assign AA clusters Done")
}else if(obj@params$method == "archetypal"){
message("LOG: main | Performing Archetypal")
obj <- obj_performArchetypal(obj, doparallel = FALSE)
message("LOG: main | Performing Archetypal Done")

message("LOG: main | assign AA clusters")
obj <- obj_assignArchetypalClusters(obj)
message("LOG: main | assign AA clusters Done")
}

# SEURAT CLUSTERIZATOIN -------------------------------
message("LOG: main | Seurat Clustering")
obj <- obj_seuratCluster(obj)
message("LOG: main | Seurat Clustering Done")

# VISUALIZE -------------------------------------------
# Visualize Dataset
message("LOG: main | Visualizing Data")
obj <- obj_visualizeData(obj)
message("LOG: main | Visualizing Data Done")

# Visualize Archetypes
if(obj@params$method=="archetypes"){
message("LOG: main | Visualizing Archetypes")
obj <- obj_visualizeArchetypes(obj)
message("LOG: main | Visualizing Archetypes Done")

# Umap Archetypes Plot
message("LOG: main | Umap Archetypes")
obj <- obj_umapArchetypes(obj)
message("LOG: main | Umap Archetypes Done")

# Archetypes Analysis
message("LOG: main | Analysis Archetypes")
obj_analysisArchetypes(obj)
message("LOG: main | Analysis Archetypes Done")
}else if(obj@params$method=="archetypal"){
message("LOG: main | Visualizing Archetypal")
obj <- obj_visualizeArchetypal(obj)
message("LOG: main | Visualizing Archetypal Done")

# Umap Archetypal Plot
message("LOG: main | Umap Archetypal")
obj <- obj_umapArchetypal(obj)
message("LOG: main | Umap Archetypal Done")

# Archetypal Analysis
message("LOG: main | Analysis Archetypal")
obj_analysisArchetypal(obj)
message("LOG: main | Analysis Archetypal Done")
}

# # Umap With Archetypes Plot
# message("LOG: main | Umap with Archetypes")
# obj <- obj_umapWithArchetypes(obj)
# message("LOG: main | Umap with Archetypes Done")
#
# obj@se@meta.data$seurat_clusters
# obj@plots$umap_orig.ident <- DimPlot(obj@se, reduction = "umap", group.by = "orig.ident")
# obj@plots$umap_tumor <- DimPlot(obj@se, reduction = "umap", group.by = "tumor")
# obj@plots$umap_seucl <- DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
# obj@plots$umap_aacl <- DimPlot(obj@se, reduction = "umap", group.by = "aa_clusters")
#
# Compare aa_clusters and seurat_clasters with Ident()
obj@compare$aa.orig <- table(obj@se@meta.data$aa_clusters, obj@se@meta.data$orig.ident)
obj@compare$se.orig <- table(obj@se@meta.data$seurat_clusters, obj@se@meta.data$orig.ident)
obj@compare$aa.se <- table(obj@se@meta.data$aa_clusters, obj@se@meta.data$seurat_clusters)

# Display this results
message("LOG: main | Comparison of clusters")
print(obj@compare$aa.orig)
print(obj@compare$se.orig)
print(obj@compare$aa.se)

# FIXME manage save obj and in particular name has been removed
# message("LOG: main | Saving Object")
# obj_saveObj(obj, name = obj@params$name)
# message("LOG: main | Saving Object Done")

