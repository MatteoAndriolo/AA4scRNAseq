# Name: allPathw
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")
# nworkers <- parallel::detectCores() - 2
# plan("multicore", workers = nworkers)

classname <- Sys.getenv("CLASSNAME")
hvf <- as.logical(Sys.getenv("HVF", "TRUE"))
test <- as.logical(Sys.getenv("TEST", "FALSE"))
test_genes <- as.numeric(Sys.getenv("TEST_GENES", 300))
test_samples <- as.numeric(Sys.getenv("TEST_SAMPLES", 500))
max_iterations <- as.numeric(Sys.getenv("MAX_ITERATIONS", 100))
num_restarts <- as.numeric(Sys.getenv("NUM_RESTARTS", 10))
out_path <- Sys.getenv("OUT_PATH")
pathw <- NULL

checkEnv <- function() {
  message("test: ", test)
  message("hvf: ", hvf)
  message("test_genes: ", test_genes)
  message("test_samples: ", test_samples)
  message("classname", classname)
  message("pathw: ", pathw)
  message("out_path: ", out_path)
  message("num_restarts: ", num_restarts)
  message("max_iterations: ", max_iterations)
}

# source("/app/Rmd/z_tools.R")
checkEnv()

obj <- new(classname)
debug <- TRUE
obj <- obj_updateParams(obj,
  updateCurrent = TRUE,
  hvf = hvf,
  test = test,
  test_genes = test_genes,
  test_samples = test_samples,
  max_iterations = max_iterations,
  num_restarts = num_restarts,
  data_path = NULL,
  out_path = out_path
)

pathw <- NULL

obj <- obj_loadData(obj, test = test, pathw = pathw, test_genes = test_genes, test_samples = test_samples)


aa.pipeline <- function(obj, pathw) {
  message("LOG: starting with PATHW ", pathw)
  obj <- obj_updateParams(obj, pathw = pathw)
  obj <- obj_loadData(obj, test = test, pathw = pathw, test_genes = test_genes, test_samples = test_samples)

  # # Visualize Dataset
  # message("Visualizing Data")
  # obj <- obj_visualizeData(obj)
  # message("Visualizing Data Done")

  # # Perform Archetypes
  # message("Performing Archetypes")
  # obj <- obj_performArchetypes(obj, max_iters = max_iterations, num_restarts = num_restarts, doparallel = FALSE)
  # message("Performing Archetypes Done")

  # # Visualize Archetypes
  # message("Visualizing Archetypes")
  # obj <- obj_visualizeArchetypes(obj, out_path)
  # message("Visualizing Archetypes Done")

  # # Umap Archetypes Plot
  # message("Umap Archetypes")
  # obj <- obj_umapArchetypes(obj, out_path)
  # message("Umap Archetypes Done")

  # message("Saving Object in ", obj@params$out_path, " ", class(obj)[[1]])
  # obj_saveObj(obj)
  # message("Saving Object Done")
}

# Register parallel backend
# numCores <- detectCores()
# cl <- makeCluster(numCores)
# registerDoParallel(cl)

# Parallel execution of the pipeline for each pathway
# foreach(pathw = pathways, .packages = c("methods")) %do% {
for (pathw in pathways) {
  message("LOG: Processing pathway: ", pathw)
  aa.pipeline(obj, pathw)
}
