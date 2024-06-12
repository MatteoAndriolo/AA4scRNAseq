#!/bin/bash

# formatt all python and R Rmd files recursively

Rscript -e "library(styler); style_dir('.', style='tidyverse', recursive=TRUE)"

black **/*.py
black **/*.ipynb