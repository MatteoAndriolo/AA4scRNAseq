#!/bin/bash
# _main.sh
#################################################
# Default parameter values
##################################################
rscriptfile="/app/Rmd/main.R"
classname="Melanoma"
hvf=FALSE
test=FALSE
verbose=FALSE
pathw=0
nworkers=10
max_iterations=100
num_restarts=10
test_genes=300
test_samples=500

##################################################
# Display usage
##################################################
echo $@
usage() {
    echo "Usage: $0 [-t test] [-H hvf] [-f rscriptfile] [-p pathw] [-g test_genes] [-s test_samples] [-w nworkers] [-h] [-v] -c classname"
    echo "  -c   Class name (required)"
    echo "  -f   Path to Rmd file (default: /app/Rmd/main.R)"
    echo "  -H   High variance filter (default: FALSE)"
    echo "  -t   Test mode (default: FALSE)"
    echo "  -g   Number of test genes (default: 300)"
    echo "  -s   Number of test samples (default: 500)"
    echo "  -p   Pathway name (default: NULL)"
    echo "  -r   Number of restarts (default: 10)"
    echo "  -i   Maximum number of iterations (default: 100)"
    echo "  -h   Show this help message and exit"
    echo "  -v   Verbose mode (default: FALSE)"
    echo "  -w   Number of workers (default: 1)"
}

##################################################
# Parse command line arguments
##################################################
while getopts ":t:g:s:H:f:p:c:r:i:v:w:h" opt; do
    case "${opt}" in
        t)
            test=${OPTARG}
            echo "Test set to: ${test}"
            ;;
        g)
            test_genes=${OPTARG}
            echo "Test genes set to: ${test_genes}"
            ;;
        s)
            test_samples=${OPTARG}
            echo "Test samples set to: ${test_samples}"
            ;;
        H)
            hvf=${OPTARG}
            echo "HVF set to: ${hvf}"
            ;;
        f)
            rscriptfile=${OPTARG}
            echo "R script file set to: ${rscriptfile}"
            ;;
        p)
            pathw="${OPTARG}"
            echo "Pathw set to: ${pathw}"
            ;;
        c)
            classname="${OPTARG}"
            echo "Class name set to: ${classname}"
            ;;
        r)
            num_restarts=${OPTARG}
            echo "Number of restarts set to: ${num_restarts}"
            ;;
        i)
            max_iterations=${OPTARG}
            echo "Max iterations set to: ${max_iterations}"
            ;;
        v) 
            verbose=true
            DEBUG="TRUE"
            ;;
        w)  
            nworkers=${OPTARG}
            echo "Number of workers set to: ${nworkers}"
            ;;
        h)
            echo "help invocked"
            usage
            exit 0
            ;;
        *)
            echo "Invalid option: -${OPTARG}" >&2
            usage
            exit 1
            ;;
    esac
done


shift $((OPTIND - 1))

##################################################
# Display parameter values for debugging
##################################################
if [ "$verbose" ]; then
    echo "LOG.sh:TEST=${test}"
    echo "LOG.sh:HVF=${hvf}"
    echo "LOG.sh:RSCRIPTFILE=${rscriptfile}"
    echo "LOG.sh:PATHW=${pathw}"
    echo "LOG.sh:CLASSNAME=${classname}"
    echo "LOG.sh:NUM_RESTARTS=${num_restarts}"
    echo "LOG.sh:MAX_ITERATIONS=${max_iterations}"
    echo "LOG.sh:TEST_GENES=${test_genes}"
    echo "LOG.sh:TEST_SAMPLES=${test_samples}"
    echo "LOG.sh:NWORKERS=${nworkers}"
    echo "LOG.sh:VERBOSE=${verbose}"
fi

##################################################
# Set output path based on classname
##################################################
case $classname in
    Exp1) output_path="out/AllonKleinLab/Experiment1" ;;
    Exp2) output_path="out/AllonKleinLab/Experiment2" ;;
    Exp3) output_path="out/AllonKleinLab/Experiment3" ;;
    Melanoma) output_path="out/Melanoma" ;;
    *) echo "Unknown classname: $classname" >&2; exit 1 ;;
esac

##################################################
# Set output folder
##################################################
prefix=""
if [ "$test" = "TRUE" ]; then
    prefix+="T"
fi
if [ "$hvf" = "TRUE" ]; then
    prefix+="H"
fi
if [ "$pathw" != 0 ]; then
    prefix+="${pathw}"
fi
timestamp=$(date +%m%d%H%M)
outpath="/app/$output_path/${timestamp}_${prefix}_$SLURM_JOB_ID"
echo "LOG: Output path is $outpath"
mkdir -p $outpath
mkdir -p $outpath/figures
mkdir -p $outpath/data

##################################################
# Initialize log file
##################################################
log_file="$outpath/stat.log"
touch $log_file
echo "Timestamp, CPU%, MEM%" > $log_file

##################################################
# Set parameters for R script 
# Set environment
##################################################
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
    "nworkers=$nworkers"
)

export CLASSNAME=$classname
export DEBUG=$DEBUG
export HVF=$hvf
export MAX_ITERATIONS=$max_iterations
export NUM_RESTARTS=$num_restarts
export NWORKERS=$nworkers
export OUT_PATH=$outpath
export PATHW="$pathw"
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