#!/bin/bash

# Define the arrays (lists)
Rmd/_run_batch.sh -v -t FALSE -m archetypal -I furthestsum -w 20 --mink 6 --maxk 7


# I=("furthestsum" "convexhull" "projected_convexhull" "partitioned_convexhull" "random")
# p=(1 2 3 4 5)
# 
# # Define the function to be executed
# your_function() {
#   paramI=$1
#   paramP=$2
#   # Add your function's implementation here
#   echo "Executing with init $paramI, pathw $paramP"
#     Rmd/_run_batch.sh -v -t FALSE -m archetypal -p $paramP -I $paramI -w 21 --mink 4 --maxk 12
# }
# 
# # Iterate over all combinations of parameters
# for paramI in "${I[@]}"; do
#   for paramP in "${p[@]}"; do
#     # Call the function with the current combination of parameters
#     your_function "$paramI" "$paramP"
#   done
# done
