library(archetypal)

## old.performArchetypal ----
# Method to perform archetypal analysis with (archetypal)
setGeneric("obj_performArchatypal", function(obj, kappas = NULL, k = NULL, doparallel = TRUE) {
  standardGeneric("obj_performArchatypal")
})

setMethod("obj_performArchatypal", "database", function(obj, kappas = NULL, k = NULL, doparallel = FALSE) {
  if (debug) {
    message("DEBUG: obj_performArchatypal | k=", k)
    message("DEBUG: obj_performArchatypal | kappas=", kappas)
    message("DEBUG: obj_performArchatypal | obj@params$k=", obj@params$k)
    message("DEBUG: obj_performArchatypal | obj@params$kappas=", obj@params$kappas)
  }
  if (is.null(obj@params$kappas) & is.null(obj@params$k)) {
    stop("ERROR: obj_performArchatypal | k and kappas are null")
  }
  if (is.null(obj@params$kappas) & !is.null(obj@params$k)) {
    obj@params$kappas <- obj@params$k
  }

  if (is.null(obj@params$which.aa)) {
    obj_updateParams(obj, which.aa = "robust")
  }

  # SETUP matrices
  message("LOG: obj_performArchatypal | Performing Archetypes on pathw ", obj@params$pathw)
  message("LOG: obj_performArchatypal | Number of archetypes is ", obj@params$kappas)

  # IMPORTANT !!!! column <-> features, row <-> samples
  m <- data.frame(t(obj_getSeData(obj)))
  colnames(m) <- rownames(obj@se)
  rownames(m) <- colnames(obj@se)

  obj@archetypes$aa.kappas <- list()

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
      message("LOG: obj_performArchatypal | Starting rerun ", i, "/", obj@params$num_restarts, " with k=", k)

      aa <- archetypal(df, opt_kappas$optimal_kappas, method = method, rseed = rseed + i * k, save_history = TRUE, nworkers = nworkers)

      history.restarts.k[[i]] <- list(
        aa = aa,
      )
    }

    t <- data.frame(
      varexpt = sapply(history.restarts.k, function(x) x$aa$varexpl),
      sse = sapply(history.restarts.k, function(x) x$aa$SSE)
    )
    t$rank_varexpl <- rank(data$varexpl, ties.method = "first")
    t$rank_sse <- rank(-data$sse, ties.method = "first")
    t$rank_sum <- data$rank_varexpl + data$rank_sse
    print(t)
    best <- which.min(data$rank_sum)
    message("INFO: obj_performArchetypal | selected ", best)

    history.restarts.k$best <- history.restarts.k[[best]]
    obj@archetypes$aa.kappas[[as.character(k)]] <- history.restarts.k
  }
  # }
  tendReruns <- Sys.time()
  message("OUTPUT: obj_performArchatypal | Reruns completed in ", difftime(tendReruns, tstartReruns, units = "secs"), " seconds")

  ## Now find the best model across all kappas
  # best_overall_run <- NULL
  # best_overall_rss <- Inf
  #
  # for (k in obj@params$kappas) {
  #  best_k_run <- obj@archetypes$aa.kappas[[as.character(k)]]$best.run
  #  if (best_k_run$rss < best_overall_rss) {
  #    best_overall_rss <- best_k_run$rss
  #    best_overall_run <- best_k_run
  #  }
  # }
  ## obj@archetypes$bestrun <- obj@archetypes$restarts[[which.min(sapply(obj@archetypes$restarts, function(x) x$rss))]]
  # obj@archetypes$bestrun <- best_overall_run
  # obj@archetypes$model <- best_overall_run$a
  # if (debug) message("DEBUG: obj_performArchatypal | dim archetypes ", dim(parameters(best_overall_run$a))[1], " ", dim(parameters(best_overall_run$a))[2])
  ## obj@archetypes$screeplot <- screeplot()
  ## obj@archetypes$model <- obj@archetypes$bestrun$a
  #
  ## obj@archetypes$restarts <- list()

  return(obj)
})
