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
TEST=FALSE
HVF=FALSE
TEST_genes=300
TEST_samples=500
CLASSNAME="Melanoma"
pathw="MAPK signaling pathway"
out_path="/app/out/Melanoma"
num_restarts=2
max_iterations=3

#docker run --rm myrocker:latest
# run unique.R script inside docker
command="Rscript /app/Rmd/unique.R"
#command="R"
docker run \
    --rm \
    -it \
    -v $(pwd):/app \
    -w /app \
    -e TEST="${TEST}" \
    -e HVF="${HVF}" \
    -e TEST_genes="${TEST_genes}" \
    -e TEST_samples="${TEST_samples}" \
    -e CLASSNAME="${CLASSNAME}" \
    -e pathw="${pathw}" \
    -e out_path="${out_path}" \
    -e num_restarts="${num_restarts}" \
    -e max_iterations="${max_iterations}" \
    -e DISABLE_AUTH=true \
    -e PASSWORD="psw" \
    myrubuntu:latest \
    $command 