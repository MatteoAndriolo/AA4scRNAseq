#!/bin/bash

# Directory containing the files to be sbatch'ed
DIRECTORY="./singularity/factory"

# Check if the directory exists
if [ ! -d "$DIRECTORY" ]; then
  echo "Directory $DIRECTORY does not exist."
  exit 1
fi

# Loop over each file in the directory
for FILE in "$DIRECTORY"/*; 
do
  # Check if it is a file
  if [ -f "$FILE" ]; then
    echo "Submitting $FILE"
    sbatch "$FILE"
  else
    echo "$FILE is not a file, skipping."
  fi
done