source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")

obj <- readRDS("/app/out/Melanoma/Melanoma_files/latestUniqueH/MelanomaH_unique.rds")
# obj@archetypes$aa.kappas=list()
# saveRDS(obj, "/app/out/Melanoma/Melanoma_files/latestUniqueH/smallHunique.rds")
# obj <- readRDS("/app/out/Melanoma/Melanoma_files/latestUniqueH/smallHunique.rds")
classname <- "Melanoma"
class(obj) <- classname
# Extract and transpose archetype parameters

#####
# Initialize a vector to store the best RSS values overall
bestrssoverall <- c()
archetype_names <- c()

# Loop through each name in the aa.kappas list
for (name in names(obj@archetypes$aa.kappas)) {
  cat("Processing name =", name, "\n")
  rss_values <- c()
  time_values <- c()

  # Loop through each element in the list corresponding to the current name
  for (j in seq_along(obj@archetypes$aa.kappas[[name]])) {
    time_value <- obj@archetypes$aa.kappas[[name]][[j]]$time
    rss_value <- obj@archetypes$aa.kappas[[name]][[j]]$rss

    cat("name =", name, ", j =", j, ", time =", time_value, ", rss =", rss_value, "\n")

    rss_values <- c(rss_values, rss_value)
    time_values <- c(time_values, time_value)
  }

  best_rss <- obj@archetypes$aa.kappas[[name]]$best.run$rss
  best_time <- obj@archetypes$aa.kappas[[name]]$best.run$time
  bestrssoverall <- c(bestrssoverall, best_rss)
  archetype_names <- c(archetype_names, name)
  cat("Best RSS of name =", name, "is", best_rss, "in time", best_time, "\n")

  mean_rss <- mean(rss_values)
  std_rss <- sd(rss_values)

  cat("Mean RSS of name =", name, "is", mean_rss, "\n")
  cat("Std Dev of RSS of name =", name, "is", std_rss, "\n\n")
}

# Find and print the overall best RSS and its corresponding index
best_rss_overall <- min(bestrssoverall)
best_rss_index <- which.min(bestrssoverall)

cat("Best RSS overall is", best_rss_overall, "in", best_rss_index, "\n")


# Create a data frame for plotting
plot_data <- data.frame(
  Archetype = archetype_names,
  Best_RSS = bestrssoverall
)

# Plot the best RSS for each archetype name using points and lines
ggplot(plot_data, aes(x = Archetype, y = Best_RSS, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Best RSS for Each Archetype",
    x = "Archetype",
    y = "Best RSS"
  ) # +
# theme(axis.text.x = element_text(angle = 45, hjust = 1))

###### with times -------
# Initialize vectors to store the best RSS values and corresponding times
bestrssoverall <- c()
archetype_names <- c()

# Initialize a list to store all RSS and time values for each archetype
all_rss_times <- list()

# Loop through each name in the aa.kappas list
for (archetype_name in names(obj@archetypes$aa.kappas)) {
  cat("Processing archetype =", archetype_name, "\n")
  rss_values <- c()
  time_values <- c()

  # Loop through each element in the list corresponding to the current archetype_name
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
ggplot(plot_data_rss, aes(x = Archetype, y = Best_RSS, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Best RSS",
    x = "Archetypes",
    y = "RSS"
  ) #+
# scale_y_continuous(expand = c(0, 0))

# Combine all time values into a single data frame for plotting
plot_data_time <- do.call(rbind, lapply(names(all_rss_times), function(name) {
  cbind(Archetype = name, all_rss_times[[name]])
}))

# Plot the time for each run grouped by archetype
ggplot(plot_data_time, aes(x = Archetype, y = Time)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Run Times",
    x = "Archetypes",
    y = "Time (sec)"
  ) +
  scale_y_continuous(limits = c(0, NA))

# Plot the RSS for each run grouped by archetype
ggplot(plot_data_time, aes(x = Archetype, y = RSS)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "RSS",
    x = "Archetypes",
    y = "RSS"
  ) #+
# scale_y_continuous(limits = c(0, NA))
