# Name: unique
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")

CLASS.NAME <- Sys.getenv("CLASSNAME")
HVF <- as.logical(Sys.getenv("HVF", "TRUE"))
TEST <- as.logical(Sys.getenv("TEST", "FALSE"))
TEST_genes <- as.numeric(Sys.getenv("TEST_genes", 300))
TEST_samples <- as.numeric(Sys.getenv("TEST_samples", 500))
max_iterations <- as.numeric(Sys.getenv("max_iterations", 100))
num_restarts <- as.numeric(Sys.getenv("num_restarts", 10))
out_path <- Sys.getenv("out_path")
pathw <- Sys.getenv("pathw")
debug <- TRUE

source("/app/Rmd/z_tools.R")
# list2env(gen_testEnv(), envir = .GlobalEnv)
# pathw <- NULL
checkEnv()

obj <- new(CLASS.NAME)
obj <- obj_updateParams(obj,
  updateCurrent = TRUE,
  HVF = HVF,
  TEST = TEST,
  TEST_genes = TEST_genes,
  TEST_samples = TEST_samples,
  max_iterations = max_iterations,
  num_restarts = num_restarts,
  data_path = NULL,
  out_path = out_path
)

obj <- obj_loadData(obj, test = TEST, pathw = pathw, test_genes = TEST_genes, test_samples = TEST_samples)


# Visualize Dataset
message("Visualizing Data")
obj <- obj_visualizeData(obj)
message("Visualizing Data Done")

# fetch number cpu
# nworkers <- parallel::detectCores()
# plan("multicore", workers = nworkers)
# Perform Archetypes
message("Performing Archetypes")
obj <- obj_performArchetypes(obj, max_iters = max_iterations, num_restarts = num_restarts, doparallel = TRUE)
message("Performing Archetypes Done")

# Visualize Archetypes
message("Visualizing Archetypes")
obj <- obj_visualizeArchetypes(obj, out_path)
message("Visualizing Archetypes Done")


# Simplex Plot
# simplexplot(obj@a)

# Umap Archetypes Plot
message("Umap Archetypes")
obj <- obj_umapArchetypes(obj, out_path)
message("Umap Archetypes Done")

message("Saving Object in ", obj@params$out_path, " ", class(obj)[[1]])
obj <- obj_saveObj(obj)
message("Saving Object Done")
