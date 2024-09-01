# Archetypal Analysis  with archetypes----

## performArchetypes ----
setGeneric("obj_performArchetypes", function(obj, kappas = NULL, k = NULL, doparallel = TRUE, doFurthestSum = TRUE) {
  standardGeneric("obj_performArchetypes")
})

setMethod("obj_performArchetypes", "database", function(obj, kappas = NULL, k = NULL, doparallel = FALSE, doFurthestSum = TRUE) {
  if (debug) {
    message("DEBUG: obj_performArchetypes | k=", k)
    message("DEBUG: obj_performArchetypes | obj@params$k=", obj@params$k)
    message("DEBUG: obj_performArchetypes | kappas=", kappas)
    message("DEBUG: obj_performArchetypes | obj@params$kappas=", obj@params$kappas)
  }
  if (is.null(obj@params$kappas) & is.null(obj@params$k)) {
    stop("ERROR: obj_performArchetypes | k and kappas are null")
  }
  if (is.null(obj@params$kappas) & !is.null(obj@params$k)) {
    obj@params$kappas <- obj@params$k
  }

  if (is.null(obj@params$which.aa)) {
    obj_updateParams(obj, which.aa = "robust")
  }

  # SETUP matrices
  message("LOG: obj_performArchetypes | Performing Archetypes on pathw ", obj@params$pathw)
  message("LOG: obj_performArchetypes | Number of archetypes is ", obj@params$kappas)

  # IMPORTANT !!!! column <-> features, row <-> samples
  m <- as.matrix(t(obj_getSeData(obj)))
  m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
  message("LOG: obj_performArchetypes | matrix dimension for archetypes is ", dim(m)[1], " ", dim(m)[2])

  # runArchetypes <- function(i, data, k, max_iterations) {
  #   message("LOG: obj_performArchetypes | Starting rerun ", i, "/", num_restarts)
  #   temp <- list()
  #   # tryCatch(
  #   #  {
  #   tstart <- Sys.time()
  #   family <- archetypes::archetypesFamily(which = "robust")
  #   temp$a <- archetypes::archetypes(data, k = k, verbose = TRUE, maxIterations = max_iterations, saveHistory = FALSE, family = family)
  #   tend <- Sys.time()
  #   message(sprintf("Archetypes Computed in %s", tend - tstart))
  #
  #   temp$rss <- temp$a$rss
  #   temp$time <- tend - tstart
  #   #  },
  #   #  error = function(e) {
  #   #    temp$a <- NULL
  #   #    temp$rss <- Inf
  #   #    temp$time <- NA
  #   #    message(sprintf("Error in archetypes: %s", e$message))
  #   #  }
  #   # )
  #   return(temp)
  # }

  # Archetypes Computation
  obj@archetypes$aa.history <- list()
  tstartReruns <- Sys.time()
  # if (doparallel) {
  #   nworkers <- parallel::detectCores() - 1
  #   # results <- mclapply(1:num_restarts, runArchetypes, data = obj@data$m, k = k, max_iterations = obj@params$max_iterations, mc.cores = 3)
  #   results <- mclapply(1:num_restarts, runArchetypes, data = m, k = k, max_iterations = obj@params$max_iterations, mc.cores = nworkers)
  #   obj@archetypes$restarts <- results
  # } else {
  family <- archetypes::archetypesFamily(which = obj@params$which.aa)
  obj@params$family <- family

  for (k in obj@params$kappas) {
    history.restarts.k <- list()
    best_rss <- Inf
    # best_model=NULL
    best_restart_index <- -1

    for (i in 1:obj@params$num_restarts) {
      temp <- list()
      tFSAstart <- Sys.time()
      message("LOG: obj_performArchetypes | Starting rerun ", i, "/", obj@params$num_restarts, " with k=", k)
      if (!is.null(doFurthestSum) & doFurthestSum) {
        obj_updateParams(obj, doFurthestSum = doFurthestSum)
        message("LOG: obj_performArchetypes | Performing Furthest Sum")
        ttstart <- Sys.time()
        irow <- sample(1:nrow(m), 1)
        irows <- FurthestSum(m, irow = irow, k = k)
        family$initfn <- archetypes:::make.fix.initfn(irows)
        ttend <- Sys.time()
        temp$FStime <- difftime(ttend, ttstart, units = "secs")
        temp$FSindices <- sort(irows)

        message("LOG: obj_performArchetypes | Furthest Sum Done in ", temp$FStime, " seconds")
        message("LOG: obj_performArchetypes | Furthest Sum Done found ", paste(temp$FSindices, collapse = ", "))
      }

      # If you want to use function instead of explicit code uncomment this
      # obj@archetypes$restarts[[i]] <- runArchetypes(i, data = obj@data$m, k = k, max_iterations = obj@params$max_iterations)

      tstart <- Sys.time()
      temp$a <- archetypes::archetypes(m, k = k, verbose = TRUE, maxIterations = obj@params$max_iterations, saveHistory = TRUE, family = family)
      tend <- Sys.time()
      tFSAend <- Sys.time()
      message("LOG: obj_performArchetypes | Archetypes Computed in ", difftime(tend, tstart, units = "secs"))

      temp$rss <- temp$a$rss
      temp$time <- difftime(tFSAend, tFSAstart, units = "secs")
      message("LOG: obj_performArchetypes | Archetypes Computed in total ", temp$time, " seconds")
      # },
      # error = function(e) {
      #   temp$a <- NULL
      #   temp$rss <- Inf
      #   temp$time <- NA
      #   message(sprintf("Error in archetypes: %s", e$message))
      # }
      # )
      history.restarts.k[[i]] <- temp

      # MISC
      message("DEBUG: obj_performArchetypes | best_rss is ", best_rss, " and temp$rss is ", temp$rss)
      if (is.na(temp$rss)) {
        temp$rss <- Inf
      }
      if (temp$rss < best_rss) {
        best_rss <- temp$rss
        # best_model = temp$a
        best_restart_index <- i
        if (debug) message("DEBUG: obj_performArchetypes | best_rss chosen is ", best_rss)
      }
    }
    history.restarts.k$best.run <- history.restarts.k[[best_restart_index]]
    obj@archetypes$aa.history[[as.character(k)]] <- history.restarts.k
    obj@archetypes$aa.bests[[as.character(k)]] <- history.restarts.k[[best]]
  }
  # }
  tendReruns <- Sys.time()
  message("OUTPUT: obj_performArchetypes | Reruns completed in ", difftime(tendReruns, tstartReruns, units = "secs"), " seconds")

  # FIXME - finding bestrun across all kappas is not the best choice, implement
  # # Now find the best model across all kappas
  # best_overall_run <- NULL
  # best_overall_rss <- Inf

  # for (k in obj@params$kappas) {
  #   best_k_run <- obj@archetypes$aa.history[[as.character(k)]]$best.run
  #   if (best_k_run$rss < best_overall_rss) {
  #     best_overall_rss <- best_k_run$rss
  #     best_overall_run <- best_k_run
  #   }
  # }
  # # obj@archetypes$bestrun <- obj@archetypes$restarts[[which.min(sapply(obj@archetypes$restarts, function(x) x$rss))]]
  # obj@archetypes$bestrun <- best_overall_run
  # obj@archetypes$model <- best_overall_run$a
  # if (debug) message("DEBUG: obj_performArchetypes | dim archetypes ", dim(parameters(best_overall_run$a))[1], " ", dim(parameters(best_overall_run$a))[2])
  # # obj@archetypes$screeplot <- screeplot()
  # # obj@archetypes$model <- obj@archetypes$bestrun$a

  # obj@archetypes$restarts <- list()

  return(obj)
})

## assignArchetypesClusters ----
setGeneric("obj_assignArchetypesClusters", function(obj) {
  standardGeneric("obj_assignArchetypesClusters")
})

setMethod("obj_assignArchetypesClusters", "database", function(obj) {
  # se <- obj@se
  # a <- obj@archetypes$model
  # k <- a$k
  message("LOG: obj_assignArchetypesClusters | creating aa_clusters metadata")
  # weights <- coef(obj@archetypes$model)
  weights <- coef(obj@archetypes$model)
  if (debug) message("DEBUG: obj_assignArchetypesClusters | dimension of weights is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  weights <- as.data.frame(weights)


  if (debug) message("LOG: obj_assignArchetypesClusters | dimension of meta.data is ", dim(obj@se@meta.data)[[1]], " ", dim(obj@se@meta.data)[[2]])
  obj@se@meta.data$aa_clusters <- apply(weights, 1, which.max)
  message("LOG: obj_assignArchetypesClusters | finished aa_clusters metadata")
  return(obj)
})

# ANALYSIS ########################
setGeneric("obj_analysisArchetypes", function(obj, ...) {
  standardGeneric("obj_analysisArchetypes")
})

setMethod("obj_analysisArchetypes", "database", function(obj, ...) {
  # Initialize vectors to store the best RSS values, times, and archetype names
  bestrssoverall <- c()
  archetype_names <- c()
  all_rss_times <- list()

  # Loop through each name in the aa.kappas list to gather RSS and time values
  for (archetype_name in names(obj@archetypes$aa.history)) {
    cat("Processing archetype =", archetype_name, "\n")
    rss_values <- c()
    time_values <- c()

    for (j in seq_along(obj@archetypes$aa.history[[archetype_name]])) {
      time_value <- obj@archetypes$aa.history[[archetype_name]][[j]]$time
      rss_value <- obj@archetypes$aa.history[[archetype_name]][[j]]$rss

      cat("archetype =", archetype_name, ", run =", j, ", time =", time_value, ", RSS =", rss_value, "\n")

      rss_values <- c(rss_values, rss_value)
      time_values <- c(time_values, time_value)
    }

    best_rss <- obj@archetypes$aa.history[[archetype_name]]$best.run$rss
    best_time <- obj@archetypes$aa.history[[archetype_name]]$best.run$time
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
  rownames(newse) <- str_replace_all(rownames(newse), "_", "-")
  # Create a temporary Seurat object with the combined matrix
  combined_obj <- CreateSeuratObject(counts = newse)
  combined_obj <- SetAssayData(combined_obj, layer = "scale.data", new.data = newse)
  # combined_obj <- ScaleData(combined_obj, layer = "counts")
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
  path_figures <- paste0(obj@params$out_path, "/figures/")

  do.call(function(...) {
    for (plot_name in names(list(...))) {
      ggsave(filename = file.path(path_figures, paste0(plot_name, ".png")), plot = obj@plots[[plot_name]])
    }
  }, results$plots)

  results

  # Save the results to an RDS file
  namefile <- paste0("resultsAnalysisunique_", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")
  output_file_path <- file.path(obj@params$out_path, namefile)
  message("LOG: obj_archetypesAnalysis | saving results in ", output_file_path)
  saveRDS(results, output_file_path)
})

## visualizeArchetypes ----
setGeneric("obj_visualizeArchetypes", function(obj) {
  standardGeneric("obj_visualizeArchetypes")
})

setMethod("obj_visualizeArchetypes", "database", function(obj) {
  # a <- obj@archetypes$model
  # k <- a$k
  path_figures <- obj@params$path_figures
  plotarchetyps <- xyplot(obj@archetypes$model, as.matrix(obj_getSeData(obj)))
  obj@plots$xyplot <- plotarchetyps
  ggsave(filename = file.path(path_figures, "archetypes_xyplot.png"), plot = obj@plots$xyplot)
  return(obj)
})

## umapArchetypes ----
setGeneric("obj_umapArchetypes", function(obj, treshold = 0.1) {
  standardGeneric("obj_umapArchetypes")
})

setMethod("obj_umapArchetypes", "database", function(obj, treshold = 0.01) {
  if (debug) message("DEBUG: obj_umapArchetypes | entering function ")
  umap_result <- UMAPPlot(obj@se)
  umap_data <- as.data.frame(umap_result$data)[, 1:2]
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  if (debug) message("DEBUG: obj_umapArchetypes | umapPlot done and fetched data")

  weights <- coef(obj@archetypes$model)
  weights <- as.data.frame(weights)
  if (debug) message("DEBUG: obj_umapArchetypes | weights dimension is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  weights[weights < treshold] <- 0

  # column_sums <- colSums(obj@archetypes$model$archetypes)
  # if (debug) message("DEBUG: obj_umapArchetypes | column_sums dimension is ", dim(column_sums)[[1]], " ", dim(column_sums)[[2]])
  # normalized_mat <- sweep(obj@archetypes$model$archetypes, 2, column_sums, FUN = "/")
  # if (debug) message("DEBUG: obj_umapArchetypes | normalized_mat dimension is ", dim(normalized_mat)[[1]], " ", dim(normalized_mat)[[2]])
  # weights <- as.data.frame(normalized_mat)

  if (debug) message("DEBUG: obj_umapArchetypes | weights dimension is ", dim(weights)[[1]], " ", dim(weights)[[2]])
  if (debug) message("DEBUG: obj_umapArchetypes | umap_data dimension is ", dim(umap_data)[[1]], " ", dim(umap_data)[[2]])

  plot_list <- list()

  for (i in 1:obj@params$k) {
    if (debug) message("DEBUG: obj_umapArchetypes | weights dimension is ", length(weights[, i]))
    umap_data$weight <- weights[, i]
    plot_title <- sprintf("UMAP Archetype %d", i)
    umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = weight)) +
      geom_point(size = 1) +
      scale_color_gradient(low = "grey", high = "red") +
      ggtitle(plot_title) +
      labs(color = "Weight")
    # print(umap_plot)

    plot_list[[i]] <- umap_plot
  }

  combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
  obj@plots$umap_archetypes <- combined_plot

  obj <- obj_updateParams(obj, updateCurrent = TRUE, umap_threshold = treshold)
  return(obj)
})

## umapWithArchetypes ----
setGeneric("obj_umapWithArchetypes", function(obj, treshold = 0.1) {
  standardGeneric("obj_umapWithArchetypes")
})

setMethod("obj_umapWithArchetypes", "database", function(obj, treshold = 0.01) {
  message("LOG: obj_umapWithArchetypes | creating plot")
  if (debug) message("DEBUG: obj_umapWithArchetypes | entering function ")
  if (debug) message("DEBUG: obj_umapWithArchetypes | archetypes dimension is ", dim(parameters(obj@archetypes$model))[[1]], " ", dim(parameters(obj@archetypes$model))[[2]])
  aspe <- t(parameters(obj@archetypes$model))
  rownames(aspe) <- rownames(obj@se@assays$RNA$counts)
  colnames(aspe) <- paste0("Archetype", 1:ncol(aspe))

  # Combine the original matrix and archetypes
  newse <- cbind(as.matrix(obj@se@assays$RNA$counts), aspe)
  newse <- as(newse, "dgCMatrix")
  if (debug) message("DEBUG: obj_umapWithArchetypes | newse dimension is ", dim(newse)[[1]], " ", dim(newse)[[2]])

  # Create a temporary Seurat object with the combined matrix
  combined_obj <- CreateSeuratObject(counts = newse)
  # TODO check if this must be done!!!
  # combined_obj <- ScaleData(combined_obj, layer = "counts")
  if (debug) message("DEBUG: obj_umapWithArchetypes | objdim ", dim(combined_obj)[[1]], " ", dim(combined_obj)[[2]])
  combined_obj <- RunPCA(combined_obj, features = rownames(combined_obj))

  # UMAP on combined matrix
  combined_obj <- RunUMAP(combined_obj, dims = 1:20) # Adjust dimensions as needed

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
  plot2 <- ggplot(combined_umap_df, aes(x = umap_1, y = umap_2, color = type)) +
    geom_point(data = subset(combined_umap_df, !type %in% archetype_labels), size = 1) +
    geom_point(data = subset(combined_umap_df, type %in% archetype_labels), color = archetype_color, size = 4) +
    geom_text(
      data = subset(combined_umap_df, type %in% archetype_labels), aes(label = as.numeric(gsub("Archetype", "", type))),
      color = "black", size = 3
    ) + # vjust = -1.5, size = 3) +
    scale_color_manual(values = c(cell_type_colors, rep(archetype_color, length(archetype_labels)))) +
    theme_minimal() +
    labs(title = "UMAP Projection of Combined SE and Archetypes", x = "UMAP 1", y = "UMAP 2")

  obj@plots$umap_with_archetypes <- plot2
  message("LOG: obj_umapWithArchetypes | finished plot")
  return(obj)
})
