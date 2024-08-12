# archetypal 2
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

    for (i in 1:obj@params$num_restarts) {
      message("LOG: obj_performArchetypal | Starting rerun ", i, "/", obj@params$num_restarts, " with k=", k)

      aa <- archetypal(df, kappas = k, method = obj@params$init_method, rseed = obj@params$rseed + i * k, save_history = FALSE, nworkers = obj@params$nworkers)

      history.restarts.k[[as.character(i)]] <- list(
        aa = aa
      )
    }

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

  t <- list()
  for (k in names(obj@archetypes$aa.bests)) {
    model <- obj@archetypes$aa.bests[[k]]
    weights <- as.data.frame(obj@archetypes$aa.bests[[k]]$A)
    # tt <- apply(weights, 1, which.max)
    tt <- apply(weights, 1, function(row) {
      max_val <- max(row)
      if (max_val > 0.6) {
        return(which.max(row))
      } else {
        return("NA")
      }
    })


    obj@archetypes$aa.bests[[k]]$cluster.id <- tt
    t[[k]] <- tt
    # obj@se@meta.data$aa_clusters <- obj@archetypes$aa.bests[[k]]$cluster.id
    obj@se$AA_clusters <- tt

    # tplot <- DimPlot(obj@se, reduction = "umap", group.by = "AA_clusters") + ggtitle("UMAP labelled by ")
    unique_tt <- unique(tt)
    num_clusters <- length(unique_tt) - 1 # excluding "NA"
    palette <- c("NA" = "grey", setNames(hue_pal()(num_clusters), unique_tt[unique_tt != "NA"]))

    tplot <- DimPlot(obj@se, reduction = "umap", group.by = "AA_clusters") +
      ggtitle("UMAP labelled by archetype if weight>0.6") +
      scale_color_manual(values = palette)

    ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", sprintf("%02d", as.numeric(k)), "_archetypes.png")), tplot)


    # PLOT ALSO ARCHETYPES POINTS
    umap_coordinates <- Embeddings(obj@se, "umap")
    archetypes_umap <- t(weights %*% umap_coordinates)

    archetypes_df <- data.frame(
      x = archetypes_umap[, 1],
      y = archetypes_umap[, 2],
      label = paste0("Archetype ", 1:nrow(archetypes_umap))
    )

    # Plot UMAP with archetypes
    tplot <- DimPlot(obj@se, reduction = "umap", group.by = "AA_clusters") +
      ggtitle("UMAP labelled by AA_clusters") +
      scale_color_manual(values = palette) +
      geom_point(data = archetypes_df, aes(x = x, y = y), color = "black", size = 3) +
      geom_text(data = archetypes_df, aes(x = x, y = y, label = label), vjust = -1, color = "red")

    ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", sprintf("%02d", as.numeric(k)), "_archLabelaAndArch.png")), tplot)
  }
  obj@se$AA_clusters <- NULL

  obj@other$AA_clusters <- t
  message("LOG: obj_assignArchetypalClusters | finished aa_clusters metadata")
  return(obj)
})

plotWithSeOrg <- function(obj, k) {
  # Get ARCHETYPES - archetypes -> names_archetypes
  archetypes <- obj@archetypes$aa.bests[[k]]$BY
  names_archetypes <- paste0("Archetype", 1:nrow(archetypes))
  rownames(archetypes) <- names_archetypes
  # Get WEIGHTS - aa.weights
  aa.weights <- obj@archetypes$aa.bests[[k]]$A
  t <- rbind(aa.weights, diag(dim(aa.weights)[2]))
  rownames(t) <- c(rownames(aa.weights), names_archetypes)
  aa.weights <- t
  aa.weights[aa.weights < treshold] <- 0

  # UMAP WITH LABELS GIVEN BY ARCHETYPES OVER TRESHOLD#########################
  obj@se.org$AA_clusters <- obj@archetypes$aa.bests[[k]]$cluster.id
  unique_aaclusters <- unique(obj@se$AA_clusters)

  # tplot <- DimPlot(obj@se, reduction = "umap", group.by = "AA_clusters") + ggtitle("UMAP labelled by ")

  num_clusters <- length(unique_aaclusters) - 1 # excluding null
  palette <- c("NA" = "grey", setNames(hue_pal()(num_clusters), unique_aaclusters[unique_tt != "NA"]))

  plotname <- paste0("UMAP_AA_", sprintf("%02d", as.numeric(k)), "_archLabels_org.png")
  tplot <- DimPlot(obj@se.org, reduction = "umap", group.by = "AA_clusters") +
    ggtitle(paste0("UMAP labelled by archetype if weight>", treshold)) +
    scale_color_manual(values = palette)

  ggsave(filename = file.path(obj@params$path_figures, plotname), tplot)

  umap_coordinates_org <- Embeddings(obj@se.org, "umap")
  archetypes_umap_org <- t(weights %*% umap_coordinates_org)

  archetypes_df_org <- data.frame(
    x = archetypes_umap_org[, 1],
    y = archetypes_umap_org[, 2],
    label = paste0("Archetype ", 1:nrow(archetypes_umap))
  )

  # Plot UMAP with archetypes
  tplot <- DimPlot(obj@se.org, reduction = "umap", group.by = "AA_clusters") +
    ggtitle("UMAP labelled by AA_clusters") +
    scale_color_manual(values = palette) +
    geom_point(data = archetypes_df_org, aes(x = x, y = y), color = "black", size = 3) +
    geom_text(data = archetypes_df_org, aes(x = x, y = y, label = label), vjust = -1, color = "red")

  ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", k, "_archLabelaAndArch_org.png")), tplot)

  # USING AL DATA
  # archetypes=as.data.frame(t(archetypes))
  tm <- GetAssayData(obj@se)
  print(tm[1:5, 1:5])
  tm <- cbind(tm, as(t(archetypes), "dgCMatrix"))

  rownames(tm) <- str_replace_all(rownames(tm), "_", "-")
  newse <- CreateSeuratObject(counts = tm)
  # newse <- NormalizeData(newse, scale.factor = 1)
  newse <- SetAssayData(object = newse, layer = "scale.data", new.data = as.matrix(tm))
  # newse <- ScaleData(newse, features = rownames(newse), do.scale = FALSE, do.center = FALSE)
  newse <- RunPCA(newse, features = rownames(newse), layers = "counts", seed.use = obj@params$rseed)
  newse <- RunUMAP(newse, features = rownames(newse), seed.use = obj@params$rseed)
  ctype <- as.vector(obj@se$ctype)
  newse$ctype <- c(ctype, rep.int(99, nrow(archetypes)))
  newse$ctype <- c(ctype, names_archetypes)

  emb <- as.data.frame(Embeddings(newse@reductions$umap))
  emb$ctypes <- factor(newse$ctype, levels = unique(newse$ctype))
  emb$rown <- rownames(newse)
  archetype_color <- "yellow"
  ctype_colors <- scales::hue_pal()(length(unique(emb$ctypes)))


  # plot <- ggplot(emb, aes(x = umap_1, y = umap_2, color = ctypes)) +
  #  geom_point(data = subset(emb, !ctypes %in% rownames(archetypes)), size = 1) +
  #  geom_point(data = subset(emb, ctypes %in% rownames(archetypes)), color = archetype_color, size = 4) +
  #  geom_text(
  #    data = subset(emb, ctypes %in% rownames(archetypes)), aes(label = as.numeric(gsub("Archetype", "", ctypes))),
  #    color = "black", size = 3
  #  ) + # vjust = -1.5, size = 3) +
  #  scale_color_manual(values = c(ctype_colors, rep(archetype_color, length(rownames(archetypes))))) +
  #  theme_minimal() +
  #  labs(title = "UMAP Projection of Combined SE and Archetypes", x = "UMAP 1", y = "UMAP 2")
  plot <- ggplot(emb, aes(x = umap_1, y = umap_2, color = ctypes)) +
    geom_point(data = subset(emb, ctypes != 99), size = 1) +
    geom_point(data = subset(emb, ctypes == 99), color = "black", size = 4) +
    geom_text(
      data = subset(emb, ctypes == 99), aes(label = as.numeric(gsub("Archetype", "", emb$rown))),
      color = "white", size = 3
    ) +
    scale_color_manual(values = c(ctype_colors, rep(archetype_color, length(names_archetypes)))) +
    theme_minimal() +
    labs(title = "CellTypes and Archetypes found", x = "UMAP 1", y = "UMAP 2")
  plot
  ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", sprintf("%02d", as.numeric(k)), ".png")), plot = plot)


  # WEIGHTS
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
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
  ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", sprintf("%02d", as.numeric(k)), "_weights.png")), plot = combined_plot)

  return(obj)
}

# obj_visualizeArchetypal -----
setGeneric("obj_visualizeArchetypal", function(obj, treshold = 0.5) {
  standardGeneric("obj_visualizeArchetypal")
})

setMethod("obj_visualizeArchetypal", "database", function(obj, treshold = 0.5) {
  message("LOG: obj_visualizeArchetypal | Visualizing archetypes")
  #  obj@archetypes$aa.history[[as.character(k)]] <- history.restarts.k
  #  obj@archetypes$aa.bests[[as.character(k)]] <- history.restarts.k[[best]]

  for (k in names(obj@archetypes$aa.bests)) {
    archetypes <- obj@archetypes$aa.bests[[k]]$BY

    names_archetypes <- paste0("Archetype", 1:nrow(archetypes))
    rownames(archetypes) <- names_archetypes

    aa.weights <- obj@archetypes$aa.bests[[k]]$A
    t <- rbind(aa.weights, diag(dim(aa.weights)[2]))
    rownames(t) <- c(rownames(aa.weights), names_archetypes)
    aa.weights <- t
    aa.weights[aa.weights < treshold] <- 0

    # FROM ASSIGN #########################
    obj@se$AA_clusters <- obj@archetypes$aa.bests[[k]]$cluster.id
    obj@se.org$AA_clusters <- obj@archetypes$aa.bests[[k]]$cluster.id

    # tplot <- DimPlot(obj@se, reduction = "umap", group.by = "AA_clusters") + ggtitle("UMAP labelled by ")
    unique_aaclusters <- unique(obj@se$AA_clusters)
    num_clusters <- length(unique_aaclusters) - 1 # excluding "NA"
    palette <- c("NA" = "grey", setNames(hue_pal()(num_clusters), unique_aaclusters[unique_tt != "NA"]))

    tplot <- DimPlot(obj@se.org, reduction = "umap", group.by = "AA_clusters") +
      ggtitle(paste0("UMAP labelled by archetype if weight>", treshold)) +
      scale_color_manual(values = palette)

    ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", k, "_archLabels_org.png")), tplot)

    tplot <- DimPlot(obj@se, reduction = "umap", group.by = "AA_clusters") +
      ggtitle(paste0("UMAP labelled by archetype if weight>", treshold)) +
      scale_color_manual(values = palette)

    ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", k, "_archLabels.png")), tplot)

    # Retrieve UMAP coordinates of the original data points
    umap_coordinates <- Embeddings(obj@se, "umap")

    # Calculate the UMAP positions for the archetypes
    archetypes_umap <- t(weights %*% umap_coordinates)

    # Create data frame for archetypes points
    archetypes_df <- data.frame(
      x = archetypes_umap[, 1],
      y = archetypes_umap[, 2],
      label = paste0("Archetype ", 1:nrow(archetypes_umap))
    )

    # Plot UMAP with archetypes
    tplot <- DimPlot(obj@se, reduction = "umap", group.by = "AA_clusters") +
      ggtitle("UMAP labelled by AA_clusters") +
      scale_color_manual(values = palette) +
      geom_point(data = archetypes_df, aes(x = x, y = y), color = "black", size = 3) +
      geom_text(data = archetypes_df, aes(x = x, y = y, label = label), vjust = -1, color = "red")

    ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", k, "_archLabelaAndArch.png")), tplot)

    #
    umap_coordinates_org <- Embeddings(obj@se.org, "umap")
    archetypes_umap_org <- t(weights %*% umap_coordinates_org)

    archetypes_df_org <- data.frame(
      x = archetypes_umap_org[, 1],
      y = archetypes_umap_org[, 2],
      label = paste0("Archetype ", 1:nrow(archetypes_umap))
    )

    # Plot UMAP with archetypes
    tplot <- DimPlot(obj@se.org, reduction = "umap", group.by = "AA_clusters") +
      ggtitle("UMAP labelled by AA_clusters") +
      scale_color_manual(values = palette) +
      geom_point(data = archetypes_df_org, aes(x = x, y = y), color = "black", size = 3) +
      geom_text(data = archetypes_df_org, aes(x = x, y = y, label = label), vjust = -1, color = "red")

    ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", k, "_archLabelaAndArch_org.png")), tplot)

    # USING AL DATA
    # archetypes=as.data.frame(t(archetypes))
    tm <- GetAssayData(obj@se)
    print(tm[1:5, 1:5])
    tm <- cbind(tm, as(t(archetypes), "dgCMatrix"))

    rownames(tm) <- str_replace_all(rownames(tm), "_", "-")
    newse <- CreateSeuratObject(counts = tm)
    # newse <- NormalizeData(newse, scale.factor = 1)
    newse <- SetAssayData(object = newse, layer = "scale.data", new.data = as.matrix(tm))
    # newse <- ScaleData(newse, features = rownames(newse), do.scale = FALSE, do.center = FALSE)
    newse <- RunPCA(newse, features = rownames(newse), layers = "counts", seed.use = obj@params$rseed)
    newse <- RunUMAP(newse, features = rownames(newse), seed.use = obj@params$rseed)
    ctype <- as.vector(obj@se$ctype)
    newse$ctype <- c(ctype, rep.int(99, nrow(archetypes)))
    newse$ctype <- c(ctype, names_archetypes)

    emb <- as.data.frame(Embeddings(newse@reductions$umap))
    emb$ctypes <- factor(newse$ctype, levels = unique(newse$ctype))
    emb$rown <- rownames(archetypes)
    archetype_color <- "yellow"
    ctype_colors <- scales::hue_pal()(length(unique(emb$ctypes)))


    # plot <- ggplot(emb, aes(x = umap_1, y = umap_2, color = ctypes)) +
    #  geom_point(data = subset(emb, !ctypes %in% rownames(archetypes)), size = 1) +
    #  geom_point(data = subset(emb, ctypes %in% rownames(archetypes)), color = archetype_color, size = 4) +
    #  geom_text(
    #    data = subset(emb, ctypes %in% rownames(archetypes)), aes(label = as.numeric(gsub("Archetype", "", ctypes))),
    #    color = "black", size = 3
    #  ) + # vjust = -1.5, size = 3) +
    #  scale_color_manual(values = c(ctype_colors, rep(archetype_color, length(rownames(archetypes))))) +
    #  theme_minimal() +
    #  labs(title = "UMAP Projection of Combined SE and Archetypes", x = "UMAP 1", y = "UMAP 2")
    plot <- ggplot(emb, aes(x = umap_1, y = umap_2, color = ctypes)) +
      geom_point(data = subset(emb, ctypes != 99), size = 1) +
      geom_point(data = subset(emb, ctypes == 99), color = "black", size = 4) +
      geom_text(
        data = subset(emb, ctypes == 99), aes(label = as.numeric(gsub("Archetype", "", emb$rown))),
        color = "white", size = 3
      ) +
      scale_color_manual(values = c(ctype_colors, rep(archetype_color, length(names_archetypes)))) +
      theme_minimal() +
      labs(title = "CellTypes and Archetypes found", x = "UMAP 1", y = "UMAP 2")
    plot
    ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", sprintf("%02d", as.numeric(k)), ".png")), plot = plot)


    # WEIGHTS
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
    combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
    ggsave(filename = file.path(obj@params$path_figures, paste0("UMAP_AA_", sprintf("%02d", as.numeric(k)), "_weights.png")), plot = combined_plot)
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
