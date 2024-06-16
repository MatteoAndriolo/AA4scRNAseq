#!/bin/bash

# Source the utility functions
if [ -f /app/Rmd/_tools.sh ]; then
    source /app/Rmd/_tools.sh
else
    source ./Rmd/_tools.sh
fi

# Default parameter values
test=TRUE
test_genes=300
test_samples=500
hvf=FALSE
RMDFILE="/app/Rmd/unique.Rmd"
pathw="MAPK signaling pathway"
classname="Melanoma"
max_iterations=100
num_restarts=10

# Parse and set parameters
set_parameters "$@"

# Validate classname and set outpath
outpath=$(set_output_path $classname)
if [ $? -ne 0 ]; then
    exit 1
fi
echo "Output path is $outpath"
mkdir -p $outpath

# Initialize log
LOG_FILE="$outpath/RMSstat.log"
initialize_log $LOG_FILE

# Render R Markdown file
#output_html="$classname.html"

# Create a parameter array to pass to render_rmd
params=(
    "rmd_file=$RMDFILE"
    "output_html=$classname.html"
    "outpath=$outpath"
    "test=$test"
    "hvf=$hvf"
    "test_genes=$test_genes"
    "test_samples=$test_samples"
    "classname=$classname"
    "pathw=$pathw"
    "max_iterations=$max_iterations"
    "num_restarts=$num_restarts"
)

render_rmd "${params[@]}"

PID=$!
echo "PID is $PID"

# Monitor the process
monitor_process $PID $LOG_FILE
