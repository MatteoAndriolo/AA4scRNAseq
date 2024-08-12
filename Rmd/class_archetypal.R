library(archetypal)

# PERFORM ARCHETYPAL ------
# Method to perform archetypal analysis with (archetypal)
setGeneric("obj_performArchetypal", function(obj, kappas = NULL, k = NULL, doparallel = TRUE) {
  standardGeneric("obj_performArchetypal")
})

setMethod("obj_performArchetypal", "database", function(obj, kappas = NULL, k = NULL, doparallel = FALSE) {
  if (debug) {
    message("DEBUG: obj_performArchetypal | k=", k)
    message("DEBUG: obj_performArchetypal | kappas=", kappas)
    message("DEBUG: obj_performArchetypal | obj@params$k=", obj@params$k)
    message("DEBUG: obj_performArchetypal | obj@params$kappas=", obj@params$kappas)
  }
  if (is.null(obj@params$kappas) & is.null(obj@params$k)) {
    stop("ERROR: obj_performArchetypal | k and kappas are null")
  }
  if (is.null(obj@params$kappas) & !is.null(obj@params$k)) {
    obj@params$kappas <- obj@params$k
  }
  # SETUP matrices
  message("LOG: obj_performArchetypal | Performing Archetypes on pathw ", obj@params$pathw)
  message("LOG: obj_performArchetypal | Number of archetypes is ", obj@params$kappas)

  # IMPORTANT !!!! column <-> features, row <-> samples
  df <- data.frame(t(GetAssayData(obj@se)))
  colnames(df) <- rownames(obj@se)
  rownames(df) <- colnames(obj@se)

  obj@archetypes$aa.bests <- list()
  obj@archetypes$aa.history <- list()

  tstartReruns <- Sys.time()
  # if (doparallel) {
  #   nworkers <- parallel::detectCores() - 1
  #   # results <- mclapply(1:num_restarts, runArchetypes, data = obj@data$m, k = k, max_iterations = obj@params$max_iterations, mc.cores = 3)
  #   results <- mclapply(1:num_restarts, runArchetypes, data = m, k = k, max_iterations = obj@params$max_iterations, mc.cores = nworkers)
  #   obj@archetypes$restarts <- results
  # } else {

  for (k in obj@params$kappas) {
    history.restarts.k <- list()
    best_rss <- Inf
    best_restart_index <- -1

    # for (i in 1:obj@params$num_restarts) {
    #   message("LOG: obj_performArchetypal | Starting rerun ", i, "/", obj@params$num_restarts, " with k=", k)

    #   aa <- archetypal(df, kappas = k, method = obj@params$init_method, rseed = obj@params$rseed + i * k, save_history = FALSE, nworkers = obj@params$nworkers)

    #   history.restarts.k[[as.character(i)]] <- list(
    #     aa = aa
    #   )
    # }
    cl2 <- makeCluster(obj@params$nworkers)
    varexport <- c("df", "k", "obj@params$init_method", "obj@params$rseed")
    clusterExport(cl2, varexport)
    history.restarts.k <- parLapply(cl2, 1:obj@params$num_restarts, function(i) {
      message("LOG: obj_performArchetypal | Starting rerun ", i, "/", obj@params$num_restarts, " with k=", k)

      aa <- archetypal(df, kappas = k, method = obj@params$init_method, rseed = obj@params$rseed + i * k, save_history = FALSE, nworkers = obj@params$nworkers)

      return(setNames(list(aa), as.character(i)))
      # history.restarts.k[[as.character(i)]] <- list(
      #   aa = aa
      # )
    })
    stopCluster(cl2)

    t <- data.frame(
      varexpt = sapply(history.restarts.k, function(x) x$aa$varexpl),
      sse = sapply(history.restarts.k, function(x) x$aa$SSE)
    )
    t$rank_varexpl <- rank(t$varexpt, ties.method = "first")
    t$rank_sse <- rank(-t$sse, ties.method = "first")
    t$rank_sum <- t$rank_varexpl + t$rank_sse
    print(t)
    best <- which.min(t$rank_sum)
    message("INFO: obj_performArchetypal | selected ", best)

    history.restarts.k$best <- history.restarts.k[[best]]$aa
    obj@archetypes$aa.history[[as.character(k)]] <- history.restarts.k
    obj@archetypes$aa.bests[[as.character(k)]] <- history.restarts.k$best
  }
  # }
  tendReruns <- Sys.time()
  message("OUTPUT: obj_performArchetypal | Reruns completed in ", difftime(tendReruns, tstartReruns, units = "secs"), " seconds")


  return(obj)
})

# obj_assignArchetypalClusters -----
setGeneric("obj_assignArchetypalClusters", function(obj) {
  standardGeneric("obj_assignArchetypalClusters")
})

setMethod("obj_assignArchetypalClusters", "database", function(obj) {
  message("LOG: obj_assignArchetypalClusters | creating aa_clusters metadata")
  #  obj@archetypes$aa.history[[as.character(k)]] <- history.restarts.k
  #  obj@archetypes$aa.bests[[as.character(k)]] <- history.restarts.k[[best]]
  treshold <- 0.5
  t <- list()
  for (k in names(obj@archetypes$aa.bests)) {
    model <- obj@archetypes$aa.bests[[k]]
    weights <- as.data.frame(obj@archetypes$aa.bests[[k]]$A)
    # tt <- apply(weights, 1, which.max)
    tt <- apply(weights, 1, function(row) {
      max_val <- max(row)
      if (max_val > treshold) {
        return(which.max(row))
      } else {
        return("NA")
      }
    })

    obj@archetypes$aa.bests[[k]]$cluster.id <- tt
    t[[k]] <- tt
    # obj@se@meta.data$aa_clusters <- obj@archetypes$aa.bests[[k]]$cluster.id
    obj@se$AA_clusters <- tt
    obj@se.org$AA_clusters <- tt

    # tplot <- DimPlot(obj@se, reduction = "umap", group.by = "AA_clusters") + ggtitle("UMAP labelled by ")
    unique_tt <- unique(tt)
    num_clusters <- length(unique_tt) - 1 # excluding "NA"
    palette <- c("NA" = "grey", setNames(hue_pal()(num_clusters), unique_tt[unique_tt != "NA"]))
    for (red in c("umap", "pca", "tsne")) {
      tplot <- DimPlot(obj@se, reduction = red, group.by = "AA_clusters") +
        ggtitle(paste0(toupper(red), " labelled by archetype if weight>0.6")) +
        scale_color_manual(values = palette)

      filename <- paste0("aa_", toupper(red), "_", sprintf("%02d", as.numeric(k)), ".png")
      ggsave(filename = file.path(obj@params$path_figures, filename), tplot)

      tplot <- DimPlot(obj@se.org, reduction = red, group.by = "AA_clusters") +
        ggtitle(paste0(toupper(red), " labelled by archetype if weight>0.5")) +
        scale_color_manual(values = palette)

      filename <- paste0("aa_", toupper(red), "_", sprintf("%02d", as.numeric(k)), "_org.png")
      ggsave(filename = file.path(obj@params$path_figures, filename), tplot)
    }
  }
  obj@se$AA_clusters <- NULL
  obj@se.org$AA_clusters <- NULL

  obj@other$AA_clusters <- t
  message("LOG: obj_assignArchetypalClusters | finished aa_clusters metadata")
  return(obj)
})

# obj_visualizeArchetypal -----
setGeneric("obj_visualizeArchetypal", function(obj, treshold = 0.5) {
  standardGeneric("obj_visualizeArchetypal")
})

setMethod("obj_visualizeArchetypal", "database", function(obj, treshold = 0.5) {
  message("LOG: obj_visualizeArchetypal | Visualizing archetypes")
  #  obj@archetypes$aa.history[[as.character(k)]] <- history.restarts.k
  #  obj@archetypes$aa.bests[[as.character(k)]] <- history.restarts.k[[best]]

  for (k in names(obj@archetypes$aa.bests)) {
    # Archetypes
    archetypes <- obj@archetypes$aa.bests[[k]]$BY
    names_archetypes <- paste0("Archetype", 1:nrow(archetypes))
    rownames(archetypes) <- names_archetypes

    # Weights
    aa.weights <- obj@archetypes$aa.bests[[k]]$A
    t <- rbind(aa.weights, diag(dim(aa.weights)[2]))
    rownames(t) <- c(rownames(aa.weights), names_archetypes)
    aa.weights <- t
    aa.weights[aa.weights < treshold] <- 0

    # Create new Seurat object with archetypes
    tm <- GetAssayData(obj@se)
    tm <- cbind(tm, as(t(archetypes), "dgCMatrix"))
    rownames(tm) <- str_replace_all(rownames(tm), "_", "-")

    newse <- CreateSeuratObject(counts = tm)
    # newse <- NormalizeData(newse, scale.factor = 1)
    # newse <- ScaleData(newse, features = rownames(newse), do.scale = FALSE, do.center = FALSE)
    newse <- SetAssayData(object = newse, layer = "scale.data", new.data = as.matrix(tm))
    newse <- RunPCA(newse, features = rownames(newse), layers = "counts", seed.use = obj@params$rseed)
    newse <- RunUMAP(newse, features = rownames(newse), seed.use = obj@params$rseed)
    newse <- RunTSNE(newse, features = rownames(newse), seed.use = obj@params$rseed, check_duplicates = FALSE)
    ctype <- as.vector(obj@se$ctype)
    newse$ctype <- c(ctype, rep.int(99, nrow(archetypes)))
    rm(ctype)

    emb <- as.data.frame(Embeddings(newse@reductions$umap))
    emb$ctype <- factor(newse$ctype, levels = unique(newse$ctype))
    emb$rown <- rownames(aa.weights)
    archetype_color <- "yellow"
    ctype_colors <- scales::hue_pal()(length(unique(emb$ctype)))

    plot <- ggplot(emb, aes(x = umap_1, y = umap_2, color = ctype)) +
      geom_point(data = subset(emb, ctype != 99), size = 1) +
      geom_point(data = subset(emb, ctype == 99), color = "black", size = 4) +
      geom_text(
        data = subset(emb, ctype == 99), aes(label = as.numeric(gsub("Archetype", "", rown))),
        color = "white", size = 3
      ) +
      scale_color_manual(values = c(ctype_colors, rep(archetype_color, length(names_archetypes)))) +
      theme_minimal() +
      labs(title = "CellTypes and Archetypes found", x = "UMAP 1", y = "UMAP 2")
    ggsave(filename = file.path(obj@params$path_figures, paste0("aa_UMAP_", sprintf("%02d", as.numeric(k)), ".png")), plot = plot)


    # Plot WEIGHTS ------
    temb <- emb
    plot_list <- list()
    for (i in 1:as.integer(k)) {
      if (debug) {
        message("DEBUG: obj_umapArchetypes | weights dimension is ", length(aa.weights[, i]))
        message("DEBUG: obj_umapArchetypes | temb dimension are ", dim(temb)[1], " ", dim(temb)[2])
      }

      temb$weight <- aa.weights[, i]
      plot_title <- sprintf("UMAP Archetype %d", i)
      umap_plot <- ggplot(temb, aes(x = umap_1, y = umap_2, color = weight)) +
        geom_point(size = 1) +
        scale_color_gradient(low = "grey", high = "red") +
        theme_minimal() +
        ggtitle(plot_title) +
        labs(color = "Weight")
      # print(umap_plot)

      plot_list[[i]] <- umap_plot
    }

    ncol <- 2
    nrow <- ceiling(as.numeric(k) / ncol)
    individual_plot_width <- 5 # width of a single plot
    individual_plot_height <- 4 # height of a single plot
    combined_plot_width <- ncol * individual_plot_width
    combined_plot_height <- nrow * individual_plot_height
    combined_plot <- plot_grid(plotlist = plot_list, ncol = ncol)
    output_filename <- file.path(obj@params$path_figures, paste0("aa_UMAP_", sprintf("%02d", as.numeric(k)), "_weights.png"))
    ggsave(filename = output_filename, plot = combined_plot, width = combined_plot_width, height = combined_plot_height)


    # HEATMAP CELL TYPE VS ARCHETYPES
    df <- data.frame(cell_type = obj@se$ctype, aa_clusters = obj@other$AA_clusters[[k]])
    contingency_table <- table(df$cell_type, df$aa_clusters)

    # Print the contingency table
    print(contingency_table)

    library(reshape2)

    # Melt the contingency table for ggplot2
    melted_table <- melt(contingency_table)

    # Plot the heatmap
    melted_plot <- ggplot(melted_table, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "blue") +
      labs(x = "Archetype", y = "Cell Type", fill = "Count") +
      theme_minimal()

    output_filename <- file.path(obj@params$path_figures, paste0("heatmap_", sprintf("%02d", as.numeric(k)), ".png"))
    ggsave(filename = output_filename, plot = melted_plot)
  }

  ## Analysis ------
  analysis_best <- data.frame(
    sse = sapply(obj@archetypes$aa.bests, function(x) x$SSE),
    varexpt = sapply(obj@archetypes$aa.bests, function(x) x$varexpl),
    time = sapply(obj@archetypes$aa.bests, function(x) x$time)
  )

  plot_best_rss <- ggplot(analysis_best, aes(x = as.numeric(rownames(analysis_best)), y = sse)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    labs(
      title = "Best RSS",
      x = "#Archetypes",
      y = "RSS"
    )
  plot_best_rss

  analysis <- list()
  for (k in names(obj@archetypes$aa.history)) {
    r <- obj@archetypes$aa.history[[k]]
    r <- r[names(r) != "best"]

    t <- data.frame(
      sse = sapply(r, function(x) x$aa$SSE),
      varexpt = sapply(r, function(x) x$aa$varexpl),
      time = sapply(r, function(x) x$aa$time)
    )
    analysis[[k]] <- t
  }
  binded <- do.call(rbind, analysis)

  plot_times <- ggplot(data = binded, aes(x = as.numeric(gsub("\\..*", "", rownames(binded))), y = time, )) +
    geom_point() +
    # geom_line() +
    theme_minimal() +
    labs(
      title = "Run Times",
      x = "#Archetypes",
      y = "Time (sec)"
    ) #+
  # scale_y_continuous(limits = c(0, NA))
  plot_times
  ggsave(filename = file.path(obj@params$path_figures, "UMAP_AA_times.png"), plot = plot_times)


  plot_sse <- ggplot(data = binded, aes(x = as.numeric(gsub("\\..*", "", rownames(binded))), y = sse, )) +
    geom_point() +
    # geom_line() +
    theme_minimal() +
    labs(
      title = "Run sse",
      x = "#Archetypes",
      y = "SSE"
    ) #+
  # scale_y_continuous(limits = c(0, NA))
  plot_sse
  ggsave(filename = file.path(obj@params$path_figures, "UMAP_AA_sse.png"), plot = plot_sse)

  plot_varexpt <- ggplot(data = binded, aes(x = as.numeric(gsub("\\..*", "", rownames(binded))), y = varexpt, )) +
    geom_point() +
    # geom_line() +
    theme_minimal() +
    labs(
      title = "Run varexpt",
      x = "#Archetypes",
      y = "varexpt"
    ) #+
  # scale_y_continuous(limits = c(0, NA))
  plot_varexpt
  ggsave(filename = file.path(obj@params$path_figures, "UMAP_AA_varexpt.png"), plot = plot_varexpt)

  # # Creating the plot with dual y-axes
  # sse_varexpt_plot <- ggplot(data = binded, aes(x = as.numeric(gsub("\\..*", "", rownames(binded)))) ) +
  #   geom_line(aes(y = sse), color = "blue") +
  #   geom_line(aes(y = varexpt), color = "red") +  # Scale varexpt for plotting
  #   scale_y_continuous(
  #     name = "SSE",
  #     sec.axis = sec_axis(~..), name = "Variance Explained")  # Scale back varexpt
  #   ) +
  #   theme_minimal() +
  #   labs(
  #     title = "Run SSE and Variance Explained",
  #     x = "#Archetypes"
  #   )
  # ggsave(filename = file.path(obj@params$path_figures, "UMAP_AA_sse_varexpt.png"), plot = sse_varexpt_plot)


  ###################################################################

  obj@archetypes$analysis <- analysis
  obj@archetypes$analysis_best <- analysis_best

  return(obj)
})
