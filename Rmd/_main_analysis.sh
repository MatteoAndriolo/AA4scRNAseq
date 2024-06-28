#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 FOLDER_PATH INPUT_FILE OUTPUT_FILE"
  exit 1
fi

# Assign command line arguments to variables
FOLDER_PATH=$1
INPUT_FILE=$2
OUTPUT_FILE=$3

# Export the environment variables for Singularity
export FOLDER_PATH=$FOLDER_PATH
export INPUT_FILE=$INPUT_FILE
export OUTPUT_FILE=$OUTPUT_FILE

# Run the R script inside the Singularity container
Rscript /app/Rmd/analysis.R
