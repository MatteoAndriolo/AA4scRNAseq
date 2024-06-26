# Name: allPathw
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")
params <- list()
params$debug <- TRUE
# params$nworkers <- parallel::detectCores() - 2
# plan("multicore", workers = nworkers)

# params$pathw <- Sys.getenv("pathw")
params$classname <- Sys.getenv("CLASSNAME")
params$hvf <- as.logical(Sys.getenv("HVF", "FALSE"))
params$max_iterations <- as.numeric(Sys.getenv("MAX_ITERATIONS", "100"))
params$num_restarts <- as.numeric(Sys.getenv("NUM_RESTARTS", "10"))
params$out_path <- Sys.getenv("OUT_PATH")
params$pathw <- as.numeric(Sys.getenv("PATHW", unset = "-1"))
params$test <- as.logical(Sys.getenv("TEST", "FALSE"))
params$test_genes <- as.numeric(Sys.getenv("TEST_GENES", "300"))
params$test_samples <- as.numeric(Sys.getenv("TEST_SAMPLES", "500"))

# TODO remove
# params$test <- TRUE
params$pathw <- -1
# params$classname <- "Melanoma"
# params$test_samples <- 2000

if (params$pathw > 0) {
  params$pathw <- pathways[[params$pathw]]
} else {
  params$pathw <- NULL
}

for (k in names(params)) {
  message("LOG: main | param ", k, ": ", params[[k]])
}
debug <- TRUE

obj <- new(params$classname)
obj <- do.call(obj_updateParams, c(list(obj = obj), params))

message("LOG: main | Loading Data")
obj <- obj_loadData(obj)
message("LOG: main | Loading Data Done")

params$pathw <- pathways[[1]]
# aa.pipeline <- function(obj, pathw) {
# for (pathw in pathways) {
message("LOG: starting with PATHW ", params$pathw)
obj <- obj_updateParams(obj, pathw = params$pathw)

message("LOG: aa.pipeline | Visualizing Data")
obj <- obj_loadData(obj)
message("LOG: aa.pipeline | Visualizing Data Done")

# Visualize Dataset
message("LOG: aa.pipeline | Visualizing Data")
obj <- obj_visualizeData(obj)
message("LOG: aa.pipeline | Visualizing Data Done")

## # fetch number cpu
## # nworkers <- parallel::detectCores()
## # plan("multicore", workers = nworkers)

# Perform Archetypes
message("LOG: aa.pipeline | Performing Archetypes")
obj <- obj_performArchetypes(obj, doparallel = FALSE)
message("LOG: aa.pipeline | Performing Archetypes Done")

# Visualize Archetypes
message("LOG: aa.pipeline | Visualizing Archetypes")
obj <- obj_visualizeArchetypes(obj)
message("LOG: aa.pipeline | Visualizing Archetypes Done")

# Umap Archetypes Plot
message("LOG: aa.pipeline | Umap Archetypes")
obj <- obj_umapArchetypes(obj)
message("LOG: aa.pipeline | Umap Archetypes Done")

message("LOG: aa.pipeline | assign AA clusters")
obj <- obj_assignAAClusters(obj)
message("LOG: aa.pipeline | assign AA clusters Done")

# message("LOG: aa.pipeline | Saving Object")
# obj_saveObj(obj)
# message("LOG: aa.pipeline | Saving Object Done")

message("LOG: aa.pipeline | LOG completed")

# Seurat Clusterizatoin
message("LOG: aa.pipeline | Seurat Clustering")
obj <- obj_seuratCluster(obj)
message("LOG: aa.pipeline | Seurat Clustering Done")
obj@se@meta.data$seurat_clusters
obj@plots$umap_orig.ident <- DimPlot(obj@se, reduction = "umap", group.by = "orig.ident")
obj@plots$umap_tumor <- DimPlot(obj@se, reduction = "umap", group.by = "tumor")
obj@plots$umap_seucl <- DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
obj@plots$umap_aacl <- DimPlot(obj@se, reduction = "umap", group.by = "aa_clusters")

# Compare aa_clusters and seurat_clasters with Ident()
obj@compare$aa.orig <- table(obj@se@meta.data$aa_clusters, obj@se@meta.data$orig.ident)
obj@compare$se.orig <- table(obj@se@meta.data$seurat_clusters, obj@se@meta.data$orig.ident)
obj@compare$aa.se <- table(obj@se@meta.data$aa_clusters, obj@se@meta.data$seurat_clusters)

message("LOG: aa.pipeline | Saving Object")
obj_saveObj(obj)
message("LOG: aa.pipeline | Saving Object Done")
# }



# Register parallel backend
# numCores <- detectCores()
# cl <- makeCluster(numCores)
# registerDoParallel(cl)

# Parallel execution of the pipeline for each pathway
# foreach(pathw = pathways, .packages = c("methods")) %do% {

# for (pathw in pathways) {
#  message("LOG: Processing pathway: ", pathw)
#  aa.pipeline(obj, pathw)
# }
