#!/bin/bash

classname="Melanoma"
RMDFILE="/app/Rmd/unique.Rmd"
test=FALSE
hvf=TRUE
test_genes=300
test_sample=500
genes=NULL
out_path=NULL

docker exec agitated_nobel Rscript -e "rmarkdown::render('$RMDFILE', params=list(TEST=$test, HVF=$hvf, TEST_genes=$test_genes, TEST_sample=$test_sample, CLASS.NAME=$classname, GENES=$genes,out_path=$outpat ))" 