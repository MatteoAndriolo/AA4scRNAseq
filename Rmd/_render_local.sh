#!/bin/bash

md=$1

docker exec agitated_nobel Rscript -e "rmarkdown::render('/app/$md')"