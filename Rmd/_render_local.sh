#!/bin/bash

# Usage: $0 [-t test] [-g test_genes] [-s test_samples] [-H hvf] [-f rscriptfile] [-p pathw] -c classname"
#   -t   Test mode (default: FALSE)"
#   -g   Number of test genes (default: 300)"
#   -s   Number of test samples (default: 500)"
#   -H   High variance filter (default: FALSE)"
#   -f   Path to Rmd file (default: /app/Rmd/unique.R)"
#   -p   Pathway name (default: NULL)"
#   -c   Class name (required)"
#   -r   Number of restarts (default: 10)"
#   -i   Maximum number of iterations (default: 100)"
#   -h   Show this help message and exit"
PARAMS=$@
echo $PARAMS


# docker exec agitated_nobel /app/Rmd/_main.sh Melanoma
command="docker run --rm -v $(pwd):/app -w /app myrubuntu /app/Rmd/_main.sh $PARAMS"
echo $command
$command