#!/bin/bash

test=TRUE
hvf=TRUE

RMDFILE="/app/Rmd/unique.Rmd"
test_genes=300
test_samples=500
genes=NULL

# Read RMD file as first parameter and extract name
classname=$1
echo $classname

case $classname in
    Exp1)
        #RMDFILE="/app/Rmd/Exp1.Rmd"
        output="out/AllonKleinLab/Experiment1"
        #classname="Exp1"
        ;;
    Exp2)
        #RMDFILE="/app/Rmd/Exp2.Rmd"
        output="out/AllonKleinLab/Experiment2"
        #classname="Exp2"
        ;;
    Exp3)
        #RMDFILE="/app/Rmd/Exp3.Rmd"
        output="out/AllonKleinLab/Experiment3"
        #classname="Exp3"
        ;;
    Melanoma)
        #RMDFILE="/app/Rmd/Melanoma.Rmd"
        output="out/Melanoma"
        #classname="Melanoma"
        ;;
    *)
        echo "Unknown classname: $classname"
        exit 1
        ;;
esac

outpath="/app/$output/${classname}_files"
output_html="$classname.html"

# create log gileo
mkdir -p $outpath
LOG_FILE="$outpath/RMSstat.log"
touch $LOG_FILE
echo "Timestamp, CPU%, MEM%" > $LOG_FILE

# -----------------------------------------------------------------
Rscript -e "rmarkdown::render('$RMDFILE', output_file='$output_html', output_dir='$outpath', params=list(TEST=$test, HVF=$hvf, TEST_genes=$test_genes, TEST_samples=$test_samples, CLASS.NAME='$classname', GENES=$genes, out_path='$outpath'))" &
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