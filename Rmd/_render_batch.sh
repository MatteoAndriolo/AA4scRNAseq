#!/bin/bash


# Create list with strings Melanoma Exp1 Exp2 Exp3 then run dfor iterating over those strings instead of using directories
listCLassNames=("Melanoma" "Exp1" "Exp2" "Exp3")

# Loop over each file in the directory
for CLASSNAME in "${listCLassNames[@]}"; do
  # Check if it is a file
  if [ -f "$CLASSNAME" ]; then
    echo "Submitting $CLASSNAME"
    sbatch "$CLASSNAME"
  else
    echo "$CLASSNAME is not a file, skipping."
  fi
done