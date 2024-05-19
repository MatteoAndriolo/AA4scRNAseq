#!/bin/bash
filepath=$1
docker run --rm -v .:/app -e DISABLE_AUTH=true -e PASSWORD=psw myubuntu/archetypes:latest Rscript -e "rmarkdown::render('$filepath')"