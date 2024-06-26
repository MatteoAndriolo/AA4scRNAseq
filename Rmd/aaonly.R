# Name: aaonly 
source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")
params <- list()
params$debug <- TRUE

params$nworkers <- parallel::detectCores() - 1
plan("multicore", workers = params$nworkers) 

# params$pathw <- Sys.getenv("pathw")
params$classname <- Sys.getenv("CLASSNAME")
params$hvf <- as.logical(Sys.getenv("HVF", "FALSE"))
params$max_iterations <- as.numeric(Sys.getenv("MAX_ITERATIONS", "100"))
params$num_restarts <- as.numeric(Sys.getenv("NUM_RESTARTS", "10"))
params$out_path <- Sys.getenv("OUT_PATH")
params$pathw <- as.numeric(Sys.getenv("PATHW", unset = "-1"))
params$test <- as.logical(Sys.getenv("TEST", "FALSE"))
params$test_genes <- as.numeric(Sys.getenv("TEST_GENES", "300"))
params$test_samples <- as.numeric(Sys.getenv("TEST_SAMPLES", "500"))

# TODO remove
# params$test <- TRUE
# params$pathw <- -1
# params$classname <- "Melanoma"
# params$test_samples <- 2000

if (params$pathw > 0) {
  params$pathw <- pathways[[params$pathw]]
} else {
  params$pathw <- NULL
}


for (k in names(params)) {
  message("LOG: main | param ", k, ": ", params[[k]])
}
debug <- TRUE

obj <- new(params$classname)
if(grepl("Exp",params$classname)){
  options(future.globals.maxSize= 5000*1024^2)
}
obj <- do.call(obj_updateParams, c(list(obj = obj), params))

message("LOG: main | Loading Data")
obj <- obj_updateParams(obj, pathw=NULL)
obj <- obj_loadData(obj)
message("LOG: main | Loading Data Done")

if(!is.null(params$pathw)){
  message("LOG: main | starting with PATHW ", params$pathw)
  obj <- obj_updateParams(obj, pathw = params$pathw)
  message("LOG: main | reloading data with pathw")
  obj <- obj_loadData(obj)
}

# Visualize Dataset
message("LOG: main | Visualizing Data")
obj <- obj_visualizeData(obj)
message("LOG: main | Visualizing Data Done")

## # fetch number cpu
## # nworkers <- parallel::detectCores()
## # plan("multicore", workers = nworkers)

# Perform Archetypes
message("LOG: main | Performing Archetypes")
main_tstart=Sys.time()
allarchetypes=list()
for(k in 4:6){
    obj <- obj_performArchetypes(obj, k=k, doparallel = FALSE)
    allarchetypes[[k]]=obj@archetypes
    message("OUT: performArchetypes | ", obj@archetypes$bestrun$time)
}
obj@other[[aa]]=allarchetypes

main_tend=Sys.time()
message("OUT: main | 10 runs of AA with k from 4 to 14 in ", difftime(main_tend, main_tstart, units="secs"), " seconds")
message("LOG: main | Performing Archetypes Done")

message("LOG: main | Saving Object")
obj_saveObj(obj, name="aaonly")
message("LOG: main | Saving Object Done")