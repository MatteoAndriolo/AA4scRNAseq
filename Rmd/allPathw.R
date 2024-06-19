# Name: allPathw
CLASS.NAME <- Sys.getenv("CLASSNAME")
HVF <- as.logical(Sys.getenv("HVF", "TRUE"))
TEST <- as.logical(Sys.getenv("TEST", "FALSE"))
TEST_genes <- as.numeric(Sys.getenv("TEST_genes", 300))
TEST_samples <- as.numeric(Sys.getenv("TEST_samples", 500))
max_iterations <- as.numeric(Sys.getenv("max_iterations", 100))
method <- Sys.getenv("method", "archetypes")
num_restarts <- as.numeric(Sys.getenv("num_restarts", 10))
out_path <- Sys.getenv("out_path")
pathw <- NULL

{
  pathways <- list(
    "Glycolysis / Gluconeogenesis",
    "MAPK signaling pathway"
  )
  pathways_lungo <- list(
    "Glycolysis / Gluconeogenesis",
    "MAPK signaling pathway",
    "mTOR signaling pathway",
    "Pathways in cancer",
    "TGF-beta signaling pathway"
  )


  checkEnv <- function() {
    message("TEST: ", TEST)
    message("HVF: ", HVF)
    message("TEST_genes: ", TEST_genes)
    message("TEST_samples: ", TEST_samples)
    message("CLASS.NAME: ", CLASS.NAME)
    message("pathw: ", pathw)
    message("out_path: ", out_path)
    message("num_restarts: ", num_restarts)
    message("max_iterations: ", max_iterations)
    message("method: ", method)
  }

  source("/app/Rmd/z_tools.R")
  checkEnv()

  if (CLASS.NAME == "") {
    stop("CLASS.NAME is NULL")
  }

  # Load required libraries and source files
  source("/app/Rmd/imports.R")
  source("/app/Rmd/classes.R")

  # Load Dataset
  obj <- new(CLASS.NAME)
  debug <- TRUE
  obj <- obj_updateParams(obj,
    updateCurrent = TRUE,
    HVF = HVF,
    TEST = TEST,
    TEST_genes = TEST_genes,
    TEST_samples = TEST_samples,
    max_iterations = max_iterations,
    num_restarts = num_restarts,
    # data_path = NA,
    out_path = out_path,
    # pathw = NULL,
    method = method
  )
  # obj@params$HVF <- HVF
  # obj@params$TEST <- TEST
  # obj@params$TEST_genes <- TEST_genes
  # obj@params$TEST_samples <- TEST_samples
  # obj@params$max_iterations <- max_iterations
  # obj@params$num_restarts <- num_restarts
  # obj@params$data_path <- NULL
  # obj@params$out_path <- out_path
  # obj@params$pathw <- NULL
  # obj@params$method <- method
  #
  # obj@curr.params <- obj@params

  obj <- obj_loadData(obj, test = TEST, pathw = pathw, test_genes = TEST_genes, test_samples = TEST_samples)
}

aa.pipeline <- function(pathw) {
  obj <- obj_updateParams(obj, pathw = pathw)
  obj <- obj_loadData(obj, test = TEST, pathw = pathw, test_genes = TEST_genes, test_samples = TEST_samples)

  # Visualize Dataset
  message("Visualizing Data")
  obj <- obj_visualizeData(obj)
  message("Visualizing Data Done")

  # Perform Archetypes
  message("Performing Archetypes")
  obj <- obj_performArchetypes(obj, max_iters = max_iterations, num_restarts = num_restarts, method = method)
  message("Performing Archetypes Done")

  # Visualize Archetypes
  message("Visualizing Archetypes")
  obj <- obj_visualizeArchetypes(obj, out_path)
  message("Visualizing Archetypes Done")

  # Umap Archetypes Plot
  message("Umap Archetypes")
  obj <- obj_umapArchetypes(obj, out_path)
  message("Umap Archetypes Done")

  message("Saving Object in ", obj@params$out_path, " ", class(obj)[[1]])
  obj_saveObj(obj)
  message("Saving Object Done")
}

# Register parallel backend
# numCores <- detectCores()
# cl <- makeCluster(numCores)
# registerDoParallel(cl)

# Parallel execution of the pipeline for each pathway
# foreach(pathw = pathways, .packages = c("methods")) %do% {
for (pathw in pathways) {
  message("Processing pathway: ", pathw)
  aa.pipeline(pathw)
}

# Stop the cluster
# stopCluster(cl)

# pathw <- pathways[[1]]
# obj <- obj_updateParams(obj, pathw = pathw)
# obj <- obj_loadData(obj, test = TEST, pathw = pathw, test_genes = TEST_genes, test_samples = TEST_samples)
#
## Visualize Dataset
# message("Visualizing Data")
# obj <- obj_visualizeData(obj)
# message("Visualizing Data Done")
#
## Perform Archetypes
# message("Performing Archetypes")
# obj <- obj_performArchetypes(obj, max_iters = max_iterations, num_restarts = num_restarts, method = method)
# message("Performing Archetypes Done")
#
## Visualize Archetypes
# message("Visualizing Archetypes")
# obj <- obj_visualizeArchetypes(obj, out_path)
# message("Visualizing Archetypes Done")
#
## Umap Archetypes Plot
# message("Umap Archetypes")
# obj <- obj_umapArchetypes(obj, out_path)
# message("Umap Archetypes Done")
#
# message("Saving Object in ", obj@params$out_path, " ", class(obj)[[1]])
# obj_saveObj(obj)
# message("Saving Object Done")
#
# pathw <- pathways[[2]]
# obj <- obj_updateParams(obj, pathw = pathw)
# obj <- obj_loadData(obj, test = TEST, pathw = pathw, test_genes = TEST_genes, test_samples = TEST_samples)
#
## Visualize Dataset
# message("Visualizing Data")
# obj <- obj_visualizeData(obj)
# message("Visualizing Data Done")
#
## Perform Archetypes
# message("Performing Archetypes")
# obj <- obj_performArchetypes(obj, max_iters = max_iterations, num_restarts = num_restarts, method = method)
# message("Performing Archetypes Done")
#
## Visualize Archetypes
# message("Visualizing Archetypes")
# obj <- obj_visualizeArchetypes(obj, out_path)
# message("Visualizing Archetypes Done")
#
## Umap Archetypes Plot
# message("Umap Archetypes")
# obj <- obj_umapArchetypes(obj, out_path)
# message("Umap Archetypes Done")
#
# message("Saving Object in ", obj@params$out_path, " ", class(obj)[[1]])
# obj_saveObj(obj)
# message("Saving Object Done")
#
