source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")

# Remove this part
debug <- TRUE
test <- FALSE
HVF <- TRUE
TEST_genes <- 300
TEST_samples <- 500
CLASS.NAME <- "Melanoma"
# pathw <- "MAPK signaling pathway"
out_path <- sprintf("/app/out/%s/%s_files", CLASS.NAME, CLASS.NAME)
obj <- new("Melanoma")
obj@params$pathw <- NULL
obj@params$test <- FALSE
obj@params$hvf <- TRUE
# read data
obj <- obj_loadData(obj)
