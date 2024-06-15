#!/bin/bash

# Parameters
test=TRUE
test_genes=300
test_samples=500
hvf=TRUE
RMDFILE="/app/Rmd/unique.Rmd"

if [ -z "$1" ]; then
    echo "No classname provided"
    exit 1
fi
classname=$1

if [ -z "$2" ]; then
    pathw=NULL
    echo "Running $classname"
else
    pathw=$2
    echo "Running $classname with pathway $pathw"
fi

# Set output path
case $classname in
    Exp1)
        output="out/AllonKleinLab/Experiment1"
        ;;
    Exp2)
        output="out/AllonKleinLab/Experiment2"
        ;;
    Exp3)
        output="out/AllonKleinLab/Experiment3"
        ;;
    Melanoma)
        output="out/Melanoma"
        ;;
    *)
        echo "Unknown classname: $classname"
        exit 1
        ;;
esac

outpath="/app/$output/${classname}_files"
mkdir -p $outpath

# Initialize log
LOG_FILE="$outpath/RMSstat.log"
touch $LOG_FILE
echo "Timestamp, CPU%, MEM%" > $LOG_FILE

# -----------------------------------------------------------------
# TEST, TEST_genes, TEST_samples, HVF, CLASS.NAME, pathw, out_path
output_html="$classname.html"

# Old - not all quoted
#Rscript -e "rmarkdown::render('$RMDFILE', output_file='$output_html', output_dir='$outpath', params=list(TEST=$test, HVF=$hvf, TEST_genes=$test_genes, TEST_samples=$test_samples, CLASS.NAME='$classname', pathw=$pathw, out_path='$outpath'))" &
# New - all quoted
Rscript -e "rmarkdown::render('$RMDFILE', output_file='$output_html', output_dir='$outpath', params=list(TEST='$test', HVF='$hvf', TEST_genes='$test_genes', TEST_samples='$test_samples', CLASS.NAME='$classname', pathw='$pathw', out_path='$outpath'))" &

PID=$!
echo "PID is $PID"

while kill -0 $PID; do
    # Fetch and log CPU and memory usage of the process
    if [ -n "$PID" ]; then
        STATS=$(top -b -n 1 | grep $PID)
        CPU=$(echo $STATS | awk '{print $5}')
        MEM=$(echo $STATS | awk '{print $9}')

        # Log to file with timestamp
        echo "$(date +%F' '%T), $CPU, $MEM" >> $LOG_FILE
    else
        echo "$(date +%F' '%T), No process found" >> $LOG_FILE
    fi
    sleep 3
done