#!/bin/bash

# Source the utility functions
if [ -f /app/Rmd/_tools.sh ]; then
    source /app/Rmd/_tools.sh
else
    source ./Rmd/_tools.sh
fi

# Default parameter values
test=FALSE
test_genes=300
test_samples=500
hvf=FALSE
RMDFILE="/app/Rmd/unique.R"
RMDFILE="/app/Rmd/allPathw.R"
pathw="NULL"
classname="Melanoma"
max_iterations=100
num_restarts=10

# Parse and set parameters
set_parameters "$@"

# Validate classname and set outpath
outpath="$(set_output_path $classname)"
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

export TEST=$test
export TEST_genes=$test_genes
export TEST_samples=$test_samples
export HVF=$hvf
export pathw=$pathw
export CLASSNAME=$classname
export max_iterations=$max_iterations
export num_restarts=$num_restarts
export out_path=$outpath


echo "classname is $classname"
#render_rmd "${params[@]}"
run_R $RMDFILE

PID=$!
echo "PID is $PID"

# Monitor the process
monitor_process $PID $LOG_FILE
