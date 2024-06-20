#!/bin/bash

# Default parameter values
classname="Melanoma"
hvf=FALSE
max_iterations=100
num_restarts=10
pathw="NULL"
rscriptfile="/app/Rmd/unique.R"
test=FALSE
test_genes=300
test_samples=500

# Display usage
usage() {
    echo "Usage: $0 [-t test] [-g test_genes] [-s test_samples] [-H hvf] [-f rscriptfile] [-p pathw] -c classname"
    echo "  -t   Test mode (default: FALSE)"
    echo "  -g   Number of test genes (default: 300)"
    echo "  -s   Number of test samples (default: 500)"
    echo "  -H   High variance filter (default: FALSE)"
    echo "  -f   Path to Rmd file (default: /app/Rmd/unique.R)"
    echo "  -p   Pathway name (default: NULL)"
    echo "  -c   Class name (required)"
    echo "  -r   Number of restarts (default: 10)"
    echo "  -i   Maximum number of iterations (default: 100)"
    echo "  -h   Show this help message and exit"
}

# Parse and set parameters
while getopts ":t:g:s:H:f:p:c:r:i:h" opt; do
    case "${opt}" in
        t) test=${OPTARG} ;;
        g) test_genes=${OPTARG} ;;
        s) test_samples=${OPTARG} ;;
        H) hvf=${OPTARG} ;;
        f) rscriptfile=${OPTARG} ;;
        p) pathw=${OPTARG} ;;
        c) classname=${OPTARG} ;;
        r) num_restarts=${OPTARG} ;;
        i) max_iterations=${OPTARG} ;;
        h)
            usage
            exit 0
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

shift $((OPTIND - 1))

# Display parameter values for debugging
echo "LOG.sh:TEST=${test}"
echo "LOG.sh:TEST_GENES=${test_genes}"
echo "LOG.sh:TEST_SAMPLES=${test_samples}"
echo "LOG.sh:HVF=${hvf}"
echo "LOG.sh:RSCRIPTFILE=${rscriptfile}"
echo "LOG.sh:PATHW=${pathw}"
echo "LOG.sh:CLASSNAME=${classname}"
echo "LOG.sh:NUM_RESTARTS=${num_restarts}"
echo "LOG.sh:MAX_ITERATIONS=${max_iterations}"

# Set output path based on classname
case $classname in
    Exp1) output_path="out/AllonKleinLab/Experiment1" ;;
    Exp2) output_path="out/AllonKleinLab/Experiment2" ;;
    Exp3) output_path="out/AllonKleinLab/Experiment3" ;;
    Melanoma) output_path="out/Melanoma" ;;
    *) echo "Unknown classname: $classname" >&2; exit 1 ;;
esac

outpath="/app/$output_path/${classname}_files"
echo "LOG: Output path is $outpath"
mkdir -p $outpath

# Initialize log file
log_file="$outpath/RMSstat.log"
touch $log_file
echo "Timestamp, CPU%, MEM%" > $log_file

params=(
    "classname=$classname"
    "hvf=$hvf"
    "max_iterations=$max_iterations"
    "num_restarts=$num_restarts"
    "outpath=$outpath"
    "output_html=$classname.html"
    "pathw=$pathw"
    "rmd_file=$rscriptfile"
    "test=$test"
    "test_genes=$test_genes"
    "test_samples=$test_samples"
)

export CLASSNAME=$classname
export HVF=$hvf
export MAX_ITERATIONS=$max_iterations
export NUM_RESTARTS=$num_restarts
export OUT_PATH=$outpath
# if pathw null do ont export

export PATHW=$pathw
export TEST=$test
export TEST_genes=$test_genes
export TEST_samples=$test_samples

# Run the R script
Rscript $rscriptfile &
pid=$!
echo "LOG.sh: PID is $pid"

# Monitor the process and log CPU and memory usage
while kill -0 $pid 2>/dev/null; do
    stats=$(top -b -n 1 | grep $pid)
    cpu=$(echo $stats | awk '{print $5}')
    mem=$(echo $stats | awk '{print $9}')
    echo "$(date +%F' '%T), $cpu, $mem" >> $log_file
    sleep 3
done
