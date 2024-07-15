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
params$method <- Sys.getenv("METHOD", "archetypal")
params$nworkers <- as.numeric(Sys.getenv("NWORKERS", "1"))
params$out_path <- Sys.getenv("OUT_PATH")
params$pathw <- as.numeric(Sys.getenv("PATHW", "0"))
params$test <- as.logical(Sys.getenv("TEST", "FALSE"))
params$test_genes <- as.numeric(Sys.getenv("TEST_GENES", "300"))
params$test_samples <- as.numeric(Sys.getenv("TEST_SAMPLES", "500"))
params$init_method <- Sys.getenv("INIT_METHOD", "furthestsum")

mink <- as.numeric(Sys.getenv("MINK"))
maxk <- as.numeric(Sys.getenv("MAXK"))
k <- as.numeric(Sys.getenv("K"))
if(!is.na(k)){
  params$kappas <- k
} else if(!(is.na(mink) | is.na(maxk))){
  params$kappas <- seq(params$mink, params$maxk)
} else{
  stop("ERROR: main | No k or mink&maxk provided")
}

debug <- params$debug
params$rseed <- 2024
if (params$pathw == 0) {
  params$pathw <- NULL
}
params$path_figures <- file.path(params$out_path, "figures")
plan("multicore", workers = params$nworkers)

################### FIXING PARAMETERS FOR TESTING
if (FALSE) {
  params$debug <- TRUE
  debug <- TRUE
  params$classname <- "Melanoma"
  params$pathw <- NULL
  params$kappas <- 2:3
  params$kappas <- 4:8
  params$test <- TRUE
  params$nworkers <- 10
  params$num_restarts <- 2
  params$max_iterations <- 5
  params$method <- "archetypal"
  params$out_path <- "/dev/null"
  params$init_method <- "furthestsum"
}
################### END FIXING PARAMETERS FOR PRESENTATION

# CREATE OBJECT -----
obj <- new(params$classname)
obj <- do.call(obj_updateParams, c(list(obj = obj), params))

# LOADING DATA -----
message("LOG: main | Loading Data")
obj <- obj_loadData(obj)
if (debug) message("DEBUG: main | Loading Data Done")

# PERFORM ARCHETYPES ---------------------------------
## ARCHETYPES----
if (obj@params$method == "archetypes") {
  message("LOG: main | Performing Archetypes")
  obj <- obj_performArchetypes(obj, doparallel = FALSE)
  message("LOG: main | Performing Archetypes Done")
 
  message("LOG: main | assign AA clusters")
  obj <- obj_assignArchetypesClusters(obj)
  message("LOG: main | assign AA clusters Done")

  ## ARCHETYPAL -----
} else if (obj@params$method == "archetypal") {
  message("LOG: main | Performing Archetypal")
  obj <- obj_performArchetypal(obj, doparallel = FALSE)
  message("LOG: main | Performing Archetypal Done")

  message("LOG: main | assign AA clusters")
  obj <- obj_assignArchetypalClusters(obj)
  message("LOG: main | assign AA clusters Done")
  obj <- obj_visualizeArchetypal(obj)
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
if (obj@params$method == "archetypes") {
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
} 


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
