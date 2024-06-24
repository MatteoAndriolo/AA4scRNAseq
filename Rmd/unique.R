# Name: unique
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
params$test <- TRUE
params$pathw <- -1
params$classname <- "Melanoma"


if (params$pathw > 0) {
  params$pathw <- pathways[[params$pathw]]
} else {
  params$pathw <- NULL
}

for (k in names(params)) {
  message("LOG: main | param ", k, ": ", params[[k]])
}

obj <- new(params$classname)

debug <- TRUE
obj <- do.call(obj_updateParams, c(list(obj = obj), params))

message("LOG: main | Loading Data")
obj <- obj_loadData(obj)
message("LOG: main | Loading Data Done")

# Visualize Dataset
message("LOG: main | Visualizing Data")
obj <- obj_visualizeData(obj)
message("LOG: main | Visualizing Data Done")

## # fetch number cpu
## # nworkers <- parallel::detectCores()
## # plan("multicore", workers = nworkers)

# Perform Archetypes
message("LOG: main | Performing Archetypes")
obj <- obj_performArchetypes(obj, doparallel = FALSE)
message("LOG: main | Performing Archetypes Done")

# Visualize Archetypes
message("LOG: main | Visualizing Archetypes")
obj <- obj_visualizeArchetypes(obj)
message("LOG: main | Visualizing Archetypes Done")

# Umap Archetypes Plot
message("LOG: main | Umap Archetypes")
obj <- obj_umapArchetypes(obj)
message("LOG: main | Umap Archetypes Done")

message("LOG: main | assign AA clusters")
obj <- obj_assignAAClusters(obj)
message("LOG: main | assign AA clusters Done")

message("LOG: main | Saving Object")
obj_saveObj(obj)
message("LOG: main | Saving Object Done")

message("LOG: main | LOG completed")

# Seurat Clusterizatoin
message("LOG: main | Seurat Clustering")
obj <- obj_seuratCluster(obj)
message("LOG: main | Seurat Clustering Done")
obj@se@meta.data$seurat_clusters
DimPlot(obj@se, reduction = "umap", group.by = "orig.ident")
DimPlot(obj@se, reduction = "umap", group.by = "tumor")
DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
DimPlot(obj@se, reduction = "umap", group.by = "aa_clusters")

# Compare aa_clusters and seurat_clasters with Ident()
table(obj@se@meta.data$aa_clusters, obj@se@meta.data$orig.ident)
table(obj@se@meta.data$seurat_clusters, obj@se@meta.data$orig.ident)
