#!/bin/bash

# Define the container name and the script to execute
# CONTAINER_NAME="agitated_nobel"
WORKDIR="/app"
SCRIPT_PATH="$WORKDIR/src/archetypes.R"

LOG_FILE="$WORKDIR/stats.log"
touch $LOG_FILE
echo "Timestamp, CPU%, MEM%" > $LOG_FILE

# Execute the Rscript inside the Docker container
Rscript $SCRIPT_PATH &
PID=$!
#PID=$(pgrep -f R) 
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


