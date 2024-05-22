#!/bin/bash

# Read RMD file as first parameter and extract name
RMDFILE=$1
echo $RMDFILE
dataset=$(basename $RMDFILE .Rmd)
echo $dataset
dataset=${dataset#singRmd}
echo $dataset

case $dataset in
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
        output="out/AllonKleinLab/Melanoma"
        ;;
    *)
        echo "Unknown filename: $dataset"
        exit 1
        ;;
esac

# create log gile
LOG_FILE="/app/$output/RMSstat.log"
touch $LOG_FILE
echo "Timestamp, CPU%, MEM%" > $LOG_FILE

# -----------------------------------------------------------------
# Run Rscript with the given parameters

Rscript -e "rmarkdown::render('$RMDFILE')" &
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


