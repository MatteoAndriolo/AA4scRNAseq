#!/bin/bash

#TEST <- as.logical(Sys.getenv("TEST", "FALSE"))
#HVF <- as.logical(Sys.getenv("HVF", "TRUE"))
#TEST_genes <- as.numeric(Sys.getenv("TEST_genes", 300))
#TEST_samples <- as.numeric(Sys.getenv("TEST_samples", 500))
#CLASS.NAME <- Sys.getenv("CLASS.NAME")
#pathw <- Sys.getenv("pathw")
#out_path <- Sys.getenv("out_path")
#num_restarts <- as.numeric(Sys.getenv("num_restarts", 10))
#max_iterations <- as.numeric(Sys.getenv("max_iterations", 100))
classname="Melanoma"
debug=TRUE
hvf=FALSE
max_iterations=3
num_restarts=2
nworkers=2
out_path="/app/out/Melanoma/$(date +%m%d_%H%M)"
mkdir -p ${out_path:5}/figures
mkdir -p ${out_path:5}/data
pathw=-1
test=TRUE
test_genes=300
test_samples=500

# run main.R script inside docker
dockercontainer=myrubuntu:latest
command="Rscript /app/Rmd/main.R"

docker run \
    --rm \
    -it \
    -v $(pwd):/app \
    -w /app \
    -e CLASSNAME="${classname}" \
    -e DEBUG="${debug}" \
    -e HVF="${hvf}" \
    -e MAX_ITERATIONS="${max_iterations}" \
    -e NUM_RESTARTS="${num_restarts}" \
    -e OUT_PATH="${out_path}" \
    -e PATHW="${pathw}" \
    -e TEST="${test}" \
    -e TEST_SAMPLES="${test_samples}" \
    -e TEST_genes="${test_genes}" \
    -e DISABLE_AUTH=true \
    -e PASSWORD="psw" \
    $dockercontainer \
    $command 