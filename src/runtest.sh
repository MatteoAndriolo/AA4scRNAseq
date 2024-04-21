#!/bin/bash

cat data/files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "File not found: {}";exit 1; fi' || exit 1 
echo "All files found"

# ensure myrocker/rstudio:archetypes exists (local docker image)
dk_arch_def="containers/Dockerfile.myrocker"
dk_arch_name="myrocker/rstudio:archetypes"

docker image inspect $dk_arch_name > /dev/null || docker build -t $dk_arch_name -f $dk_arch_def .


