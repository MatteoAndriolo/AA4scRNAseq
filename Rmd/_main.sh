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
method="archetypal"
nworkers=10
max_iterations=100
num_restarts=10
test_genes=300
test_samples=500
init_method="furthestsum"
k=0
mink=0
maxk=0

##################################################
# Display usage
##################################################
##################################################
# Display usage
##################################################
echo $@
usage() {
    echo "Usage: $0 [-t test] [-H hvf] [-f rscriptfile] [-p pathw] [-g test_genes] [-s test_samples] [-m method] [-I initmethod] [-w nworkers] [-h] [-v] -c classname"
    echo "  -c   Class name (required)"
    echo "  -f   Path to Rmd file (default: /app/Rmd/main.R)"
    echo "  -H   High variance filter (default: FALSE)"
    echo "  -t   Test mode (default: FALSE)"
    echo "  -g   Number of test genes (default: 300)"
    echo "  -s   Number of test samples (default: 500)"
    echo "  -p   Pathway name (default: NULL)"
    echo "  -r   Number of restarts (default: 10)"
    echo "  -i   Maximum number of iterations (default: 100)"
    echo "  -I   Select init method (default: furthestsum)"
    echo "  -h   Show this help message and exit"
    echo "  -v   Verbose mode (default: FALSE)"
    echo "  -m   Archetypal Analysis method (default: archetypal)"
    echo "  -w   Number of workers (default: 1)"
    echo "  -k   Number of archetypes (default: NULL)"
    echo "  --mink Minimum number of archetypes (default: NULL)"
    echo "  --maxk Maximum number of archetypes (default: NULL)"
}

##################################################
# Parse command line arguments
##################################################
PARSED_OPTIONS=$(getopt -o t:g:s:H:f:p:c:r:k:i:I:vw:m:h --long mink:,maxk: -- "$@")
if [ $? -ne 0 ]; then
    usage
    exit 1
fi

eval set -- "$PARSED_OPTIONS"

while true; do
    case "$1" in
        -t)
            test=$2
            echo "Test set to: ${test}"
            shift 2
            ;;
        -g)
            test_genes=$2
            echo "Test genes set to: ${test_genes}"
            shift 2
            ;;
        -s)
            test_samples=$2
            echo "Test samples set to: ${test_samples}"
            shift 2
            ;;
        -H)
            hvf=$2
            echo "HVF set to: ${hvf}"
            shift 2
            ;;
        -f)
            rscriptfile=$2
            echo "R script file set to: ${rscriptfile}"
            shift 2
            ;;
        -p)
            pathw=$2
            echo "Pathw set to: ${pathw}"
            shift 2
            ;;
        -c)
            classname=$2
            echo "Class name set to: ${classname}"
            shift 2
            ;;
        -r)
            num_restarts=$2
            echo "Number of restarts set to: ${num_restarts}"
            shift 2
            ;;
        -k)
            k=$2
            echo "Number of archetypes set to: ${k}"
            shift 2
            ;;
        --mink)
            mink=$2
            echo "Minimum number of archetypes set to: ${mink}"
            shift 2
            ;;
        --maxk)
            maxk=$2
            echo "Maximum number of archetypes set to: ${maxk}"
            shift 2
            ;;
        -i)
            max_iterations=$2
            echo "Max iterations set to: ${max_iterations}"
            shift 2
            ;;
        -I)
            init_method=$2
            echo "Init method set to: ${init_method}"
            shift 2
            ;;
        -v) 
            verbose=true
            DEBUG="TRUE"
            shift 1
            ;;
        -w)  
            nworkers=$2
            echo "Number of workers set to: ${nworkers}"
            shift 2
            ;;
        -m) 
            method=$2
            echo "Method selected is ${method}"
            shift 2
            ;;
        -h)
            echo "Help invoked"
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Invalid option: -$1" >&2
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
    echo "LOG: _main.sh | env | TEST=${test}"
    echo "LOG: _main.sh | env | HVF=${hvf}"
    echo "LOG: _main.sh | env | RSCRIPTFILE=${rscriptfile}"
    echo "LOG: _main.sh | env | PATHW=${pathw}"
    echo "LOG: _main.sh | env | CLASSNAME=${classname}"
    echo "LOG: _main.sh | env | NUM_RESTARTS=${num_restarts}"
    echo "LOG: _main.sh | env | MAX_ITERATIONS=${max_iterations}"
    echo "LOG: _main.sh | env | TEST_GENES=${test_genes}"
    echo "LOG: _main.sh | env | TEST_SAMPLES=${test_samples}"
    echo "LOG: _main.sh | env | METHOD=${method}"
    echo "LOG: _main.sh | env | NWORKERS=${nworkers}"
    echo "LOG: _main.sh | env | VERBOSE=${verbose}"
    echo "LOG: _main.sh | env | K=${k}"
    echo "LOG: _main.sh | env | MINK=${mink}"
    echo "LOG: _main.sh | env | MAXK=${maxk}"
fi

##################################################
# Set output path based on classname
##################################################
case $classname in
    Exp1) output_path="out/AllonKleinLab/Experiment1" ;;
    Exp2) output_path="out/AllonKleinLab/Experiment2" ;;
    Exp3) output_path="out/AllonKleinLab/Experiment3" ;;
    Melanoma) output_path="out/Melanoma" ;;
    Mouse) output_path="out/MouseCortex" ;;
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
if [ "$init_method" == "furthestsum" ]; then
    prefix+="FS"
elif [ "$init_method" == "convexhull" ]; then
    prefix+="CH"
elif [ "$init_method" == "projected_convexhull"]; then
    prefix+="PRCH"
elif [ "$init_method" == "partitioned_convexhull"]; then
    prefix+="PACH"
elif [ "$init_method" == "random"]; then
    prefix+="RND"
else

    prefix+="${init_method:0:2}"
fi
if [ "$pathw" != 0 ]; then
    prefix+="${pathw}"
fi

timestamp=$(date +%m%d_%H%M)
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
    "method=$method"
    "init_method=$init_method"
    "nworkers=$nworkers"
)

export CLASSNAME=$classname
export DEBUG=$DEBUG
export HVF=$hvf
export MAX_ITERATIONS=$max_iterations
export NUM_RESTARTS=$num_restarts
export METHOD=$method
export INIT_METHOD=$init_method
export NWORKERS=$nworkers
export OUT_PATH=$outpath
export PATHW="$pathw"
export TEST=$test
export TEST_genes=$test_genes
export TEST_samples=$test_samples
export K=$k
export MINK=$mink
export MAXK=$maxk

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