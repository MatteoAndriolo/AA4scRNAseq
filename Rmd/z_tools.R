checkEnv <- function() {
  message("TEST: ", TEST)
  message("HVF: ", HVF)
  message("TEST_genes: ", TEST_genes)
  message("TEST_samples: ", TEST_samples)
  message("CLASS.NAME: ", CLASS.NAME)
  message("pathw: ", pathw)
  message("out_path: ", out_path)
  message("num_restarts: ", num_restarts)
  message("max_iterations: ", max_iterations)
}

gen_testEnv <- function() {
  return(list(
    "TEST" = TRUE,
    "HVF" = TRUE,
    "TEST_genes" = 300,
    "TEST_samples" = 400,
    "CLASS.NAME" = "Melanoma",
    "pathw" = NULL,
    "out_path" = "/app/out",
    "num_restarts" = 3,
    "max_iterations" = 5
  ))
}

