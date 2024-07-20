#!/bin/bash

#Rmd/_run_batch.sh -t TRUE -p 1 -m archetypal -I furthestsum -w 10 --mink 4 --maxk 5 -c Melanoma  -v
#Rmd/_run_batch.sh -t TRUE -p 1 -m archetypal -I convexhull -w 10 --mink 4 --maxk 5 -c Melanoma  -v

# exit 1

# Define the arrays (lists)
# I=("furthestsum" "convexhull" "projected_convexhull" "partitioned_convexhull" "random")
# p=(1 2 3 4 5)
I=("furthestsum") 
p=(1)

# Define the function to be executed
your_function() {
  paramI=$1
  paramP=$2
  # Add your function's implementation here
  echo "Executing with init $paramI, pathw $paramP"
    Rmd/_run_batch.sh -v -t TRUE -m archetypal -p $paramP -I $paramI -w 20 --mink 4 --maxk 6
}

# Iterate over all combinations of parameters
for paramI in "${I[@]}"; do
  for paramP in "${p[@]}"; do
    # Call the function with the current combination of parameters
    your_function "$paramI" "$paramP"
  done
done
