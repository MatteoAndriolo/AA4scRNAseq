# Name: unique
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")
debug <- TRUE
params <- list()
params$classname <- Sys.getenv("CLASSNAME")
params$debug <- debug
params$hvf <- as.logical(Sys.getenv("HVF", "FALSE"))
params$k <- as.numeric(Sys.getenv("K", "8"))
params$max_iterations <- as.numeric(Sys.getenv("MAX_ITERATIONS", "100"))
params$num_restarts <- as.numeric(Sys.getenv("NUM_RESTARTS", "10"))
params$nworkers <- as.numeric(Sys.getenv("NWORKERS", "10"))
params$out_path <- Sys.getenv("OUT_PATH")
params$pathw <- as.numeric(Sys.getenv("PATHW", unset = "-1"))
params$test <- as.logical(Sys.getenv("TEST", "FALSE"))
params$test_genes <- as.numeric(Sys.getenv("TEST_GENES", "300"))
params$test_samples <- as.numeric(Sys.getenv("TEST_SAMPLES", "500"))

# params$nworkers <- parallel::detectCores() - 2
plan("multicore", workers = params$nworkers)

################### FIXING PARAMETERS FOR PRESENTATION
# TODO remove this
# FOR NOW WE WILL NOT USE PATHW BUT ONLY HVF!!!
params$pathw <- -1
if (params$pathw > 0) {
  params$pathw <- pathways[[params$pathw]]
} else {
  params$pathw <- NULL
}

# if(is.null(params$k) & classname=="Melanoma"){
params$kappas <- 7:15
# }

params$hvf <- TRUE
################### END FIXING PARAMETERS FOR PRESENTATION

for (k in names(params)) {
  message("LOG: main | param ", k, ": ", params[[k]])
}

# CREATE OBJECT -----
obj <- new(params$classname)
obj <- do.call(obj_updateParams, c(list(obj = obj), params))

# LOADING DATA -----
message("LOG: main | Loading Data")
obj <- obj_loadData(obj)
message("LOG: main | Loading Data Done")

# if(!is.null(params$pathw)){
#   message("LOG: main | starting with PATHW ", params$pathw)
#   obj <- obj_updateParams(obj, pathw = params$pathw)
#   message("LOG: main | reloading data with pathw")
#   obj <- obj_loadData(obj)
# }

# PERFORM ARCHETYPES ---------------------------------
message("LOG: main | Performing Archetypes")
obj <- obj_performArchetypes(obj, doparallel = FALSE)
message("LOG: main | Performing Archetypes Done")

message("LOG: main | assign AA clusters")
obj <- obj_assignAAClusters(obj)
message("LOG: main | assign AA clusters Done")

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
message("LOG: main | Visualizing Archetypes")
obj <- obj_visualizeArchetypes(obj)
message("LOG: main | Visualizing Archetypes Done")

# Umap Archetypes Plot
message("LOG: main | Umap Archetypes")
obj <- obj_umapArchetypes(obj)
message("LOG: main | Umap Archetypes Done")

# # Umap With Archetypes Plot
# message("LOG: main | Umap with Archetypes")
# obj <- obj_umapWithArchetypes(obj)
# message("LOG: main | Umap with Archetypes Done")
#
# obj@se@meta.data$seurat_clusters
# obj@plots$umap_orig.ident <- DimPlot(obj@se, reduction = "umap", group.by = "orig.ident")
# obj@plots$umap_tumor <- DimPlot(obj@se, reduction = "umap", group.by = "tumor")
# obj@plots$umap_seucl <- DimPlot(obj@se, reduction = "umap", group.by = "seurat_clusters")
# obj@plots$umap_aacl <- DimPlot(obj@se, reduction = "umap", group.by = "aa_clusters")
#
# # Compare aa_clusters and seurat_clasters with Ident()
# obj@compare$aa.orig <- table(obj@se@meta.data$aa_clusters, obj@se@meta.data$orig.ident)
# obj@compare$se.orig <- table(obj@se@meta.data$seurat_clusters, obj@se@meta.data$orig.ident)
# obj@compare$aa.se <- table(obj@se@meta.data$aa_clusters, obj@se@meta.data$seurat_clusters)

message("LOG: main | Saving Object")
obj_saveObj(obj, name = "unique")
message("LOG: main | Saving Object Done")

# FROM ANALYSIS --------------------------------------------------------
# Fetch folder path, input file name, and output file name from environment variables

# Check if environment variables are provided
if (params$out_path == "" || input_file == "" || output_file == "") {
  stop("params$out_path, INPUT_FILE, and OUTPUT_FILE environment variables must be provided.")
}

# Construct full paths for input and output files
input_file_path <- file.path(params$out_path, input_file)
# create namegile resultsAnalysisunique with date
namefile <- paste0("resultsAnalysisunique_", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")
output_file_path <- file.path(params$out_path, namefile)

# Initialize vectors to store the best RSS values, times, and archetype names
bestrssoverall <- c()
archetype_names <- c()
all_rss_times <- list()

# Loop through each name in the aa.kappas list to gather RSS and time values
for (archetype_name in names(obj@archetypes$aa.kappas)) {
  cat("Processing archetype =", archetype_name, "\n")
  rss_values <- c()
  time_values <- c()

  for (j in seq_along(obj@archetypes$aa.kappas[[archetype_name]])) {
    time_value <- obj@archetypes$aa.kappas[[archetype_name]][[j]]$time
    rss_value <- obj@archetypes$aa.kappas[[archetype_name]][[j]]$rss

    cat("archetype =", archetype_name, ", run =", j, ", time =", time_value, ", RSS =", rss_value, "\n")

    rss_values <- c(rss_values, rss_value)
    time_values <- c(time_values, time_value)
  }

  best_rss <- obj@archetypes$aa.kappas[[archetype_name]]$best.run$rss
  best_time <- obj@archetypes$aa.kappas[[archetype_name]]$best.run$time
  bestrssoverall <- c(bestrssoverall, best_rss)
  archetype_names <- c(archetype_names, archetype_name)

  all_rss_times[[archetype_name]] <- data.frame(Run = seq_along(rss_values), Time = time_values, RSS = rss_values)

  cat("Best RSS for archetype =", archetype_name, "is", best_rss, "in time", best_time, "\n")

  mean_rss <- mean(rss_values)
  std_rss <- sd(rss_values)

  cat("Mean RSS for archetype =", archetype_name, "is", mean_rss, "\n")
  cat("Standard Deviation of RSS for archetype =", archetype_name, "is", std_rss, "\n\n")
}

# Find and print the overall best RSS and its corresponding index
best_rss_overall <- min(bestrssoverall)
best_rss_index <- which.min(bestrssoverall)

cat("Best RSS overall is", best_rss_overall, "for archetype", archetype_names[best_rss_index], "\n")

# Create a data frame for plotting best RSS values
# Convert archetype names to numeric values
numeric_archetype_names <- as.numeric(gsub("Archetype", "", archetype_names))
plot_data_rss <- data.frame(
  Archetype = factor(numeric_archetype_names, levels = sort(numeric_archetype_names)),
  Best_RSS = bestrssoverall
)

# Plot the best RSS for each archetype name using points and lines
plot1 <- ggplot(plot_data_rss, aes(x = Archetype, y = Best_RSS, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Best RSS",
    x = "#Archetypes",
    y = "RSS"
  )

# Combine all time values into a single data frame for plotting
plot_data_time <- do.call(rbind, lapply(names(all_rss_times), function(name) {
  cbind(Archetype = as.numeric(gsub("Archetype", "", name)), all_rss_times[[name]])
}))

# Plot the time for each run grouped by archetype
plot2 <- ggplot(plot_data_time, aes(x = factor(Archetype, levels = sort(unique(Archetype))), y = Time)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Run Times",
    x = "#Archetypes",
    y = "Time (sec)"
  ) +
  scale_y_continuous(limits = c(0, NA))

# Plot the RSS for each run grouped by archetype
plot3 <- ggplot(plot_data_time, aes(x = factor(Archetype, levels = sort(unique(Archetype))), y = RSS)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "RSS",
    x = "#Archetypes",
    y = "RSS"
  )

# Extract and transpose archetype parameters
aspe <- t(parameters(obj@archetypes$model))
rownames(aspe) <- rownames(obj@se@assays$RNA$counts)
colnames(aspe) <- paste0("Archetype", 1:ncol(aspe))

# Combine the original matrix and archetypes
newse <- cbind(as.matrix(obj@se@assays$RNA$counts), aspe)
newse <- as(newse, "dgCMatrix")

# Create a temporary Seurat object with the combined matrix
combined_obj <- CreateSeuratObject(counts = newse)
combined_obj <- ScaleData(combined_obj, layer = "counts")
combined_obj <- RunPCA(combined_obj, features = rownames(combined_obj))

# UMAP on combined matrix
combined_obj <- RunUMAP(combined_obj, dims = 1:20)

# Extract UMAP embeddings
umap_combined <- Embeddings(combined_obj, "umap")

# Separate the combined results for plotting
combined_umap_df <- data.frame(umap_combined)

# Add cell types for the original cells
cell_types <- obj@se@meta.data$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK.
archetype_labels <- colnames(aspe)
combined_umap_df$type <- c(cell_types, archetype_labels)

# Convert cell types to factors to ensure consistent ordering
combined_umap_df$type <- factor(combined_umap_df$type, levels = c(unique(cell_types), archetype_labels))

# Define color for archetypes and other points
archetype_color <- "yellow"
cell_type_colors <- scales::hue_pal()(length(unique(cell_types)))

# Plot UMAP with special points highlighted
plot4 <- ggplot(combined_umap_df, aes(x = umap_1, y = umap_2, color = type)) +
  geom_point(data = subset(combined_umap_df, !type %in% archetype_labels), size = 1) +
  geom_point(data = subset(combined_umap_df, type %in% archetype_labels), color = archetype_color, size = 4) +
  geom_text(
    data = subset(combined_umap_df, type %in% archetype_labels), aes(label = as.numeric(gsub("Archetype", "", type))),
    color = "black", size = 3
  ) +
  scale_color_manual(values = c(cell_type_colors, rep(archetype_color, length(archetype_labels)))) +
  theme_minimal() +
  labs(title = "UMAP plot with Archetypes", x = "UMAP 1", y = "UMAP 2")

# Create a summary of time and RSS by number of archetypes
summary_data <- data.frame(
  Archetype = archetype_names,
  Best_RSS = bestrssoverall,
  Mean_RSS = sapply(all_rss_times, function(x) mean(x$RSS)),
  Std_Dev_RSS = sapply(all_rss_times, function(x) sd(x$RSS)),
  Mean_Time = sapply(all_rss_times, function(x) mean(x$Time)),
  Std_Dev_Time = sapply(all_rss_times, function(x) sd(x$Time)),
  Original_RSS = I(all_rss_times),
  Original_Times = I(all_rss_times)
)

# Create a list to store the plots, new data, and summary data
results <- list(
  plots = list(
    best_rss = plot1,
    run_times = plot2,
    rss_per_run = plot3,
    umap_with_archetypes = plot4
  ),
  new_data = list(
    combined_matrix = newse,
    umap_embeddings = umap_combined
  ),
  summary = summary_data
)

results

# Save the results to an RDS file
message("LOG: saving results in ", output_file_path)
saveRDS(results, output_file_path)
