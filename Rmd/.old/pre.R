# Name: unique
sink(stdout(), type = "message")
sink(stdout(), type = "output")

# Import necessary R scripts
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")

# Load required libraries
library(parallel)
library(Seurat)
library(Rtsne)
library(plotly)
library(ggplot2)
library(cowplot)
devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
install.packages("viridis")
library(viridis)

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
params$out_path <- "out/Melanoma/full"
params$test_genes <- 300
params$test_samples <- 500
params$init_method <- "furthestsum"
mink <- 6
maxk <- 18
k <- mink:maxk

mink <- 5
maxk <- 6
k <- mink:maxk

params$rseed <- 2024
set.seed(params$rseed)
params$path_figures <- file.path(params$out_path, "figures")
dir.create(params$path_figures)
params$path_outdata <- file.path(params$out_path, "data")
dir.create(params$path_outdata)

params$pathw <- NULL
params$hvf <- FALSE
params$test <- FALSE
params$test <- TRUE

obj <- new(params$classname)
obj <- do.call(obj_updateParams, c(list(obj = obj), params))
obj <- obj_loadData(obj)
obj_visualizeData(obj)


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
  "TGF-beta signaling pathway",
  NULL
)

name_pathways <- list(
  "Glycolysis / Gluconeogenesis" = "GLYC",
  "MAPK signaling pathway" = "MAPK",
  "mTOR signaling pathway" = "MTOR",
  "Pathways in cancer" = "CANC",
  "TGF-beta signaling pathway" = "TGFB",
  "HVF" = "HVF"
)

# Helper function to find closest points
findClosestPoints <- function(se, aa, k) {
  archetype <- t(aa$aa.bests[[k]]$BY)
  aa_clusters <- aa$aa.bests[[k]]$cluster.id
  aa_clusters_treshold <- aa$aa.bests[[k]]$cluster.id.treshold.5
  num_archetypes <- ncol(archetype)
  weights <- as.data.frame(aa$aa.bests[[k]]$A)
  common <- intersect(rownames(se), rownames(archetype))

  small_se <- se[common, ]

  closestPoints <- sapply(1:num_archetypes, function(i) {
    distances <- apply(GetAssayData(small_se), 2, function(cell) sum((cell[common] - archetype[, i])^2))
    which.min(distances)
  })

  newCtype <- as.character(se$ctype)
  newCtype[closestPoints] <- "Archetype"
  newCtype <- factor(newCtype, levels = c(levels(se$ctype), "Archetype"))
  weights[closestPoints, ] <- diag(num_archetypes)
  rownames(weights) <- colnames(se)
  colnames(weights) <- paste("Archetype", 1:num_archetypes, sep = "_")

  se <- FindVariableFeatures(se)

  if (!"pca" %in% Reductions(se)) {
    se <- RunPCA(se, features = VariableFeatures(se))
  }
  se <- RunUMAP(se, features = VariableFeatures(se), dim.embed = 3)
  se <- RunTSNE(se, dims = 1:30, perplexity = 30, check_duplicates = FALSE, dim.embed = 3)

  se$ctype <- newCtype
  se@misc$aa_clusters <- aa_clusters
  se@misc$aa_cluster_treshold.5 <- aa$aa.bests[[k]]$cluster.id.treshold.5
  se@misc$weights <- weights

  return(se)
}

# Create cluster
# cl <- makeCluster(params$nworkers)
cl <- makeCluster(2)
clusterExport(cl, c("params", "list_parallel_params", "pathways", "name_pathways", "findClosestPoints", "obj", "obj_updateParams", "obj_loadData", "obj_visualizeData", "obj_performArchetypal", "obj_assignArchetypalClusters", "obj_visualizeArchetypal", "obj_seuratCluster"))

# Parallel execution
# res <- future.apply::future_lapply(1:length(list_parallel_params), function(i) {
res <- parLapply(cl, 1:length(list_parallel_params), function(i) {
  tryCatch(
    {
      paramsT <- params
      paramsT$hvf <- list_parallel_params[[i]]$HVF
      paramsT$pathw <- list_parallel_params[[i]]$pathw
      if (paramsT$hvf) {
        name <- "HVF"
        paramsT$pathw <- NULL
      } else if (paramsT$pathw > 0) {
        name <- name_pathways[[pathways[[paramsT$pathw]]]]
      } else {
        stop("Unexpected values")
      }

      paramsT$out_path <- file.path(paramsT$out_path, name)
      if (!dir.exists(paramsT$out_path)) {
        dir.create(paramsT$out_path, recursive = TRUE)
      }

      obj <- do.call(obj_updateParams, c(list(obj = obj), paramsT))
      obj <- obj_loadData(obj)
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

      saveRDS(obj@other, file = file.path(paramsT$out_path, "metadata.Rds"))
      saveRDS(obj@archetypes, file = file.path(paramsT$out_path, "archetypes.Rds"))

      return(list(HVF = list_parallel_params[[i]]$HVF, pathw = list_parallel_params[[i]]$pathw, name = name, out = paramsT$out_path))
    },
    error = function(e) {
      message <- paste0("Error in iteration ", i, ": ", e$message)
      cat(message, "\n")
      return(message)
    }
  )
})

# Print results
print(res)

# Stop cluster
stopCluster(cl)
