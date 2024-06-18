#!/bin/bash

# formatt all python and R Rmd files recursively

Rscript -e "library(styler); style_dir('./Rmd', recursive=TRUE)"
Rscript -e "library(styler); style_dir('./src', recursive=TRUE)"
Rscript -e "library(styler); style_dir('./archetypes', recursive=TRUE)"

black **/*.py
black **/*.ipynb