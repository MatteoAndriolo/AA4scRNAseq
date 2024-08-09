# Name: unique
sink(stdout(), type = "message")
sink(stdout(), type = "output")
# Import necessary R scripts
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")

# Fetch params from environment variables
params <- list()
params$classname <- "Melanoma"
params$debug <- TRUE
debug <- TRUE
params$num_restarts <- "10"
params$method <- "archetypal"
params$nworkers <- 20
if (params$nworkers > parallel::detectCores()) {
  params$nworkers <- parallel::detectCores()
}
params$out_path <- "out/Melanoma"
params$test_genes <- 300
params$test_samples <- 500
params$init_method <- "furthestsum"
mink <- 6
maxk <- 18
k <- mink:maxk

params$rseed <- 2024
params$path_figures <- file.path(params$out_path, "figures")
params$path_outdata <- file.path(params$out_path, "data")

params$pathw <- NULL
params$HVF <- FALSE
params$test <- TRUE
obj <- new(params$classname)
obj <- do.call(obj_updateParams, c(list(obj = obj), params))
obj <- obj_loadData(obj)
obj_visualizeData(obj)

plan("multicore", workers = params$nworkers)
# set memory limit to 30GB
options("future.global.maxSize" = 100 * 1024^3)

list_parallel_params <- list(
  list(HVF = TRUE, pathw = NULL),
  list(HVF = FALSE, pathw = 1),
  list(HVF = FALSE, pathw = 2),
  list(HVF = FALSE, pathw = 3),
  list(HVF = FALSE, pathw = 4),
  list(HVF = FALSE, pathw = 5)
)

pathways <- list(
  "Glycolysis / Gluconeogenesis",
  "MAPK signaling pathway",
  "mTOR signaling pathway",
  "Pathways in cancer",
  "TGF-beta signaling pathway"
)

name_pathways <- list(
  "Glycolysis / Gluconeogenesis" = "GLYC",
  "MAPK signaling pathway" = "MAPK",
  "mTOR signaling pathway" = "MTOR",
  "Pathways in cancer" = "CANC",
  "TGF-beta signaling pathway" = "TGFB"
)

# ececute following code in parallel using obj_updateParams(list_parallel_params[[i]]) to update params for each different run
cl <- makeCluster(params$nworkers)
clusterExport(cl, c("params", "list_parallel_params", "obj", "pathways", "name_pathways"))


res <- parSapply(cl, 1:length(list_parallel_params), function(i) {
  paramsT <- params
  paramsT$HVF <- list_parallel_params[[i]]$HVF
  paramsT$pathw <- list_parallel_params[[i]]$pathw

  if (paramsT$HVF) {
    name <- "HVF"
  } else if (paramsT$pathw) {
    name <- name_pathways[[pathways[[paramsT$pathw]]]]
  } else {
    # return error -> values not expected
    stop("Unexpected values")
  }

  paramsT$out_path <- file.path(paramsT$out_path, name)
  if (!dir.exists(paramsT$out_path)) {
    dir.create(paramsT$out_path, recursive = TRUE)
  }

  obj <- do.call(obj_updateParams, c(list(obj = obj), paramsT))
  obj <- obj_loadData(obj) # pathw or HVF
  obj_visualizeData(obj, name = paste0("_", name))

  obj <- obj_performArchetypal(obj, doparallel = TRUE)

  obj <- obj_assignArchetypalClusters(obj)

  obj <- obj_visualizeArchetypal(obj)

  obj <- obj_seuratCluster(obj)

  obj@archetypes$seurat_clusters <- obj@se$seurat_clusters
  obj@archetypes$ctypes <- obj@se$ctype
  obj@other$se.metadata <- obj@se@meta.data
  obj@other$genenames <- rownames(obj@se)
  obj@other$cellnames <- colnames(obj@se)
  saveRDS(obj@other, file = file.path(obj@params$out_path, "metadata.Rds"))
  saveRDS(obj@arcohetypes, file = file.path(obj@params$out_path, "archetypes.Rds"))
  if (!list_parallel_params[[i]]$HVF) {
    Sys.sleep(1)
  }
  return(list(HVF = list_parallel_params[[i]]$HVF, pathw = list_parallel_params[[i]]$pathw, name = name, out = paramsT$out_path))
})

res

stopCluster(cl)
