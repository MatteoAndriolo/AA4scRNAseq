#!/bin/bash

# Directory containing the files to be sbatch'ed
DIRECTORY="./Rmd/factory"

# Check if the directory exists
if [ ! -d "$DIRECTORY" ]; then
  echo "Directory $DIRECTORY does not exist."
  exit 1
fi

# Usage: $0 [-t test] [-H hvf] [-f rscriptfile] [-p pathw] [-g test_genes] [-s test_samples] [-w nworkers] [-h] [-v] -c classname
#   -c   Class name (required)
#   -f   Path to Rmd file (default: /app/Rmd/main.R)
#   -H   High variance filter (default: FALSE)
#   -t   Test mode (default: FALSE)
#   -g   Number of test genes (default: 300)
#   -s   Number of test samples (default: 500)
#   -p   Pathway name (default: NULL)
#   -r   Number of restarts (default: 10)
#   -i   Maximum number of iterations (default: 100)
#   -h   Show this help message and exit
#   -v   Verbose mode (default: FALSE)
#   -w   Number of workers (default: 1)
PARAMS="-c Melanoma -f /app/Rmd/main.R -w 20 -v -p 1"

# Loop over each file in the directory
for FILE in "$DIRECTORY"/*; 
do
  # Check if it is a file
  if [ -f "$FILE" ]; then
    echo "Submitting $FILE"
    sbatch "$FILE" $PARAMS
  else
    echo "$FILE is not a file, skipping."
  fi
done