#!/bin/bash

# Check if the notebook file argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <notebook-file>"
  exit 1
fi

# Set the notebook file environment variable
export NOTEBOOK_FILE=$1

# Run the nbconvert command
conda run -n pyaa \
    jupyter nbconvert --to html --execute --config nbconvert_execute.py --output "${NOTEBOOK_FILE%.ipynb}.html" "$NOTEBOOK_FILE" 
