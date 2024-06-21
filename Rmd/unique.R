# Name: unique
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")
debug <- TRUE
# nworkers <- parallel::detectCores() - 2
# plan("multicore", workers = nworkers)

# pathw <- Sys.getenv("pathw")
classname <- Sys.getenv("CLASSNAME")
hvf <- as.logical(Sys.getenv("HVF", "FALSE"))
max_iterations <- as.numeric(Sys.getenv("MAX_ITERATIONS", "100"))
num_restarts <- as.numeric(Sys.getenv("NUM_RESTARTS", "10"))
out_path <- Sys.getenv("OUT_PATH")
pathw <- as.numeric(Sys.getenv("PATHW", unset = "-1"))
test <- as.logical(Sys.getenv("TEST", "FALSE"))
test_genes <- as.numeric(Sys.getenv("TEST_GENES", "300"))
test_samples <- as.numeric(Sys.getenv("TEST_SAMPLES", "500"))


if (pathw > 0) {
  pathw <- pathways[[pathw]]
} else {
  pathw <- NULL
}

message("LOG: TEST: ", test)
message("LOG: HVF: ", hvf)
message("LOG: TEST_genes: ", test_genes)
message("LOG: TEST_samples: ", test_samples)
message("LOG: CLASSNAME: ", classname)
message("LOG: PATHW: ", pathw)
message("LOG: OUT_PATH: ", out_path)
message("LOG: NUM_RESTARTS: ", num_restarts)
message("LOG: MAX_ITERATIONS: ", max_iterations)

obj <- new(classname)

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

obj <- obj_loadData(obj, test = test, pathw = pathw, test_genes = test_genes, test_samples = test_samples)

# Visualize Dataset
message("Visualizing Data")
obj <- obj_visualizeData(obj)
message("Visualizing Data Done")

## # fetch number cpu
## # nworkers <- parallel::detectCores()
## # plan("multicore", workers = nworkers)

# Perform Archetypes
message("Performing Archetypes")
obj <- obj_performArchetypes(obj, doparallel = FALSE)
message("Performing Archetypes Done")

# Visualize Archetypes
message("Visualizing Archetypes")
obj <- obj_visualizeArchetypes(obj)
message("Visualizing Archetypes Done")

# Umap Archetypes Plot
message("Umap Archetypes")
obj <- obj_umapArchetypes(obj)
message("Umap Archetypes Done")

message("Saving Object")
obj_saveObj(obj)
message("Saving Object Done")

message("LOG completed")
