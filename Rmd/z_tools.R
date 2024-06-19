gen_testEnv <- function() {
  return(list(
    "TEST" = TRUE,
    "HVF" = TRUE,
    "TEST_genes" = 300,
    "TEST_samples" = 2000,
    "CLASS.NAME" = "Melanoma",
    "pathw" = NULL,
    "out_path" = "/app/out",
    "num_restarts" = 3,
    "max_iterations" = 5,
    "method" = "archetypes"
  ))
}

# add list to env
list2env(gen_testEnv(), envir = .GlobalEnv)
