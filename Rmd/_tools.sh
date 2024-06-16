#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 [-t test] [-g test_genes] [-s test_samples] [-H hvf] [-f RMDFILE] [-p pathw] -c classname"
    echo "  -t   Test mode (default: FALSE)"
    echo "  -g   Number of test genes (default: 300)"
    echo "  -s   Number of test samples (default: 500)"
    echo "  -H   High variance filter (default: FALSE)"
    echo "  -f   Path to Rmd file (default: /app/Rmd/unique.Rmd)"
    echo "  -p   Pathway name (default: MAPK signaling pathway)"
    echo "  -c   Class name (required)"
    echo "  -r   Number of restarts (default: 10)"
    echo "  -i   Maximum number of iterations (default: 100)"
    echo "  -h   Show this help message and exit"
}

# Function to parse and set parameters
set_parameters() {
    while getopts ":t:g:s:H:f:p:c:r:i:h" opt; do
        case "${opt}" in
            t) test=${OPTARG} ;;
            g) test_genes=${OPTARG} ;;
            s) test_samples=${OPTARG} ;;
            H) hvf=${OPTARG} ;;
            f) RMDFILE=${OPTARG} ;;
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

    # Validate classname
    if [ -z "$classname" ]; then
        echo "No classname provided"
        usage
        exit 1
    fi

    # Display the values (for debugging purposes)
    echo "test=${test}"
    echo "test_genes=${test_genes}"
    echo "test_samples=${test_samples}"
    echo "hvf=${hvf}"
    echo "RMDFILE=${RMDFILE}"
    echo "pathw=${pathw}"
    echo "classname=${classname}"
    echo "num_restarts=${num_restarts}"
    echo "max_iterations=${max_iterations}"
}

# Function to set output path based on classname
set_output_path() {
    local classname=$1
    local output_path

    case $classname in
        Exp1)
            output_path="out/AllonKleinLab/Experiment1"
            ;;
        Exp2)
            output_path="out/AllonKleinLab/Experiment2"
            ;;
        Exp3)
            output_path="out/AllonKleinLab/Experiment3"
            ;;
        Melanoma)
            output_path="out/Melanoma"
            ;;
        *)
            echo "Unknown classname: $classname" >&2
            return 1
            ;;
    esac

    echo "/app/$output_path/${classname}_files"
}

# Function to initialize log file
initialize_log() {
    local log_file=$1
    touch $log_file
    echo "Timestamp, CPU%, MEM%" > $log_file
}

# Function to render the R Markdown file
render_rmd() {
    declare -A params
    for param in "$@"; do
        IFS="=" read -r key value <<< "$param"
        params[$key]=$value
    done

    Rscript -e "rmarkdown::render('${params[rmd_file]}', output_file='${params[output_html]}', output_dir='${params[outpath]}', params=list(TEST='${params[test]}', HVF='${params[hvf]}', TEST_genes='${params[test_genes]}', TEST_samples='${params[test_samples]}', CLASS.NAME='${params[classname]}', pathw='${params[pathw]}', out_path='${params[outpath]}', max_iterations='${params[max_iterations]}', num_restarts='${params[num_restarts]}'))" &
}

# Function to monitor the process and log CPU and memory usage
monitor_process() {
    local pid=$1
    local log_file=$2

    while kill -0 $pid 2>/dev/null; do
        if [ -n "$pid" ]; then
            local stats=$(top -b -n 1 | grep $pid)
            local cpu=$(echo $stats | awk '{print $5}')
            local mem=$(echo $stats | awk '{print $9}')

            # Log to file with timestamp
            echo "$(date +%F' '%T), $cpu, $mem" >> $log_file
        else
            echo "$(date +%F' '%T), No process found" >> $log_file
        fi
        sleep 3
    done
}
