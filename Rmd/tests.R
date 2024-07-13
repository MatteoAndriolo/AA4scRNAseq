source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")

# Remove this part
debug <- TRUE
CLASS.NAME <- "Melanoma"
# pathw <- "MAPK signaling pathway"
out_path <- sprintf("/app/out/%s/%s_files", CLASS.NAME, CLASS.NAME)
obj <- new("Melanoma")
obj@params$pathw <- NULL
obj@params$test <- TRUE
obj@params$hvf <- TRUE
obj@params$test_genes <- 300
obj@params$test_samples <- 500
# read data
obj <- obj_loadData(obj)

kappas <- NULL
k <- NULL
doparallel <- FALSE
obj@params$nworkers <- 20

# IMPORTANT !!!! column <-> features, row <-> samples
df <- data.frame(t(obj_getSeData(obj)))
colnames(df) <- rownames(obj@se)
rownames(df) <- colnames(obj@se)

# Find Optimal number of k -----
# runs for each kappas from 1 to max and stores SSE VARIANCE EXPLAINED from
#  each run and computer elbow with UIK
find_optimal_kappas()

# Initial approx for AA ------
# convexhull
# # find_outmost_convexhull_points
# # find_outmost_projected_convexhull_points - ch for all possible combinations of variables taken by npr (default 2)
# # find_outmost_partitioned_convexhull_points - makes np partitions of data frame then computes ch for each partition and gives ch of overall union
# outmost_points - compute outmost outermost points
# furthesum
# random


optk <- find_optimal_kappas(df, nworkers = 10)

obj@archetypes$aa.kappas <- list()
aa <- archetypal(df = df, kappas = 3, verbose = TRUE, rseed = 9102, save_history = TRUE)
aa$BY

print(aa)
summary(aa)
names(aa$run_results)
aa$run_results$time
aa$run_results
