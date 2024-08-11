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
library(ggsankey)
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
params$out_path <- "out/Melanoma"
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
params$path_outdata <- file.path(params$out_path, "data")

params$pathw <- NULL
params$HVF <- FALSE
params$test <- FALSE
params$test <- TRUE

# Initialize object
obj <- new(params$classname)
obj <- do.call(obj_updateParams, c(list(obj = obj), params))
obj <- obj_loadData(obj)
obj_visualizeData(obj)

plan("multicore", workers = params$nworkers)
plan("multisession", workers = params$nworkers)
# set memory limit to 30GB
options("future.global.maxSize" = 100 * 1024 ^ 3)

list_parallel_params <- list(
#  list(HVF = TRUE, pathw = NULL),
#  list(HVF = FALSE, pathw = 1),
#  list(HVF = FALSE, pathw = 2),
#  list(HVF = FALSE, pathw = 3),
#  list(HVF = FALSE, pathw = 4),
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

  small_se <- se[common,]

  closestPoints <- sapply(1:num_archetypes, function(i) {
    distances <- apply(GetAssayData(small_se), 2, function(cell) sum((cell[common] - archetype[, i]) ^ 2))
    which.min(distances)
  })

  newCtype <- as.character(se$ctype)
  newCtype[closestPoints] <- "Archetype"
  newCtype <- factor(newCtype, levels = c(levels(se$ctype), "Archetype"))
  weights[closestPoints,] <- diag(num_archetypes)
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
# clusterExport(cl, c("params", "list_parallel_params", "obj", "pathways", "name_pathways", "findClosestPoints"))

# Parallel execution
res <- future_apply(cl, 1:length(list_parallel_params), function(i) {
  paramsT <- params
  paramsT$HVF <- list_parallel_params[[i]]$HVF
  paramsT$pathw <- list_parallel_params[[i]]$pathw

  if (paramsT$HVF) {
    name <- "HVF"
  } else if (paramsT$pathw) {
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

  # TODO kappas for analysis
  kappas <- list("7", "12")
  kappas <- list(mink) #TODO for testing only, remove this line
  #ik <- 2
  #k <- kappas[[ik]]
  for (k in kappas) {
    num_archetypes <- as.integer(k)
    newse <- findClosestPoints(obj@se, obj@archetypes, k)

    # Setup parameters for plotting
    colors_types <- viridis(length(unique(newse$ctype)) - 1)
    size_types <- 1
    colors_Archetype <- "black"
    size_Archetype <- 4
    colors_text_Archetype <- c("white")
    size_text_Archetype <- 3

    # Setup data for plotting
    plot_data <- as.data.frame(Embeddings(newse, reduction = "tsne"))
    colnames(plot_data) <- c("X1", "X2", "X3")
    plot_data$Label <- newse$ctype
    plot_data$aa_clusters <- newse@misc$aa_clusters
    plot_data$aa_clusters_treshold.5 <- newse@misc$aa_cluster_treshold.5

    # Generate and save plots
    pairs <- list(c("X1", "X2"), c("X1", "X3"), c("X2", "X3"))
    plots <- list()

    for (p in pairs) {
      x <- p[1]
      y <- p[2]

      p2d <- ggplot(plot_data, aes_string(x = x, y = y, color = "Label")) +
      geom_point(data = subset(plot_data, Label != "Archetype"), size = size_types) +
      geom_point(data = subset(plot_data, Label == "Archetype"), color = "black", size = size_Archetype) +
      geom_text(
        data = subset(plot_data, Label == "Archetype"),
        aes(label = seq(1, num_archetypes)),
        color = colors_text_Archetype, size = size_text_Archetype
      ) +
      scale_color_manual(values = c(colors_types, colors_Archetype)) +
      theme_minimal() +
      ggtitle(paste(name, x, "vs", y))

      plots[[length(plots) + 1]] <- p2d

      ggsave(
      file.path(paramsT$out_path, paste0(name, ".tsne.", x, "vs", y, ".png")),
      p2d,
      width = 8,
      height = 6
    )
    }

    combined_plot <- plot_grid(plotlist = plots, ncol = 1)
    ggsave(
    file.path(paramsT$out_path, paste0(name, ".tsne.combined.png")),
    combined_plot,
    width = 8,
    height = 18
  )
  }
  # Save results
  saveRDS(obj@other, file = file.path(paramsT$out_path, "metadata.Rds"))
  saveRDS(obj@archetypes, file = file.path(paramsT$out_path, "archetypes.Rds"))

  return(list(HVF = list_parallel_params[[i]]$HVF, pathw = list_parallel_params[[i]]$pathw, name = name, out = paramsT$out_path))
})

# Print results
print(res)

# Stop cluster
stopCluster(cl)
