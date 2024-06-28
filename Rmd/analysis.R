# Load required libraries and source files
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")

# Load the object
obj <- readRDS("/app/out/Melanoma/Melanoma_files/latestUniqueH/MelanomaH_unique.rds")
classname <- "Melanoma"
class(obj) <- classname

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
plot_data_rss <- data.frame(
  Archetype = factor(archetype_names, levels = archetype_names),
  Best_RSS = bestrssoverall
)

# Plot the best RSS for each archetype name using points and lines
plot1 <- ggplot(plot_data_rss, aes(x = Archetype, y = Best_RSS, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Best RSS for Each Archetype",
    x = "Archetype",
    y = "Best RSS"
  )

# Combine all time values into a single data frame for plotting
plot_data_time <- do.call(rbind, lapply(names(all_rss_times), function(name) {
  cbind(Archetype = name, all_rss_times[[name]])
}))

# Plot the time for each run grouped by archetype
plot2 <- ggplot(plot_data_time, aes(x = Archetype, y = Time)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Run Times for Each Archetype",
    x = "Archetype",
    y = "Time (sec)"
  ) +
  scale_y_continuous(limits = c(0, NA))

# Plot the RSS for each run grouped by archetype
plot3 <- ggplot(plot_data_time, aes(x = Archetype, y = RSS)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "RSS for Each Run",
    x = "Archetype",
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
  labs(title = "UMAP Projection with Archetypes", x = "UMAP 1", y = "UMAP 2")

# Create a summary of time and RSS by number of archetypes
summary_data <- data.frame(
  Archetype = archetype_names,
  Best_RSS = bestrssoverall,
  Mean_RSS = sapply(all_rss_times, function(x) mean(x$RSS)),
  Std_Dev_RSS = sapply(all_rss_times, function(x) sd(x$RSS)),
  Mean_Time = sapply(all_rss_times, function(x) mean(x$Time)),
  Std_Dev_Time = sapply(all_rss_times, function(x) sd(x$Time))
)

# Add original RSS and time data to the summary
# summary_data$Original_RSS <- all_rss_times
# summary_data$Original_Times <- all_rss_times

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
  summary = summary_data,
  table_runs = all_rss_times
)

# Save the results to an RDS file
# saveRDS(results, "/app/out/Melanoma/Melanoma_files/latestUniqueH/MelanomaH_results.rds")
