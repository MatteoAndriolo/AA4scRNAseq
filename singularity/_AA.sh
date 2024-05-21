#!/bin/bash

ANALYSIS=$1
DATASET=$2
NUM_ARCHETYPES=$3

if [ -z NUM_ARCHETYPES ]; then 
    NUM_ARCHETYPES=0
fi

WORKDIR="/app"
SCRIPT_PATH="$WORKDIR/src/main.R"

LOG_FILE="$WORKDIR/out/$DATASET/stats.log"
touch $LOG_FILE
echo "Timestamp, CPU%, MEM%" > $LOG_FILE

# -----------------------------------------------------------------
# Run Rscript with the given parameters

Rscript $SCRIPT_PATH $ANALYSIS $DATASET $NUM_ARCHETYPES &
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


