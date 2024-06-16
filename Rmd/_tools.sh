#!/bin/bash

# Default values
test=FALSE
test_genes=300
test_samples=500
hvf=FALSE
RMDFILE="/app/Rmd/unique.Rmd"
pathw="MAPK signaling pathway"

# Function to display usage
usage() {
  echo "Usage: $0 [-t test] [-g test_genes] [-s test_samples] [-H hvf] [-r RMDFILE] [-p pathw]"
  echo "  -t   Test mode (default: FALSE)"
  echo "  -g   Number of test genes (default: 300)"
  echo "  -s   Number of test samples (default: 500)"
  echo "  -H   High variance filter (default: FALSE)"
  echo "  -r   Path to Rmd file (default: /app/Rmd/unique.Rmd)"
  echo "  -p   Pathway name (default: MAPK signaling pathway)"
  echo "  -h   Show this help message and exit"
}

# Function to parse and set parameters
set_parameters() {
  while getopts ":t:g:s:H:r:p:h" opt; do
    case "${opt}" in
      t)
        test=${OPTARG}
        ;;
      g)
        test_genes=${OPTARG}
        ;;
      s)
        test_samples=${OPTARG}
        ;;
      H)
        hvf=${OPTARG}
        ;;
      r)
        RMDFILE=${OPTARG}
        ;;
      p)
        pathw=${OPTARG}
        ;;
      h)
        usage
        exit 0
        ;;
      *)
        usage
        exit 1
        ;;
    esac
  done

  shift $((OPTIND-1))

  # Display the values (for debugging purposes)
  echo "test=${test}"
  echo "test_genes=${test_genes}"
  echo "test_samples=${test_samples}"
  echo "hvf=${hvf}"
  echo "RMDFILE=${RMDFILE}"
  echo "pathw=${pathw}"
}

