#!/bin/bash

# Directory containing the files to be sbatch'ed
DIRECTORY="./Rmd/factory"

# Check if the directory exists
if [ ! -d "$DIRECTORY" ]; then
  echo "Directory $DIRECTORY does not exist."
  exit 1
fi

PARAMS="-f /app/Rmd/main.R -w 20 -v TRUE -p 1  -c Melanoma -m archetypes"

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