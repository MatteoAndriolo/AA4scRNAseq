parse_parameters() {
    if [ "$#" -ne 6 ]; then
        echo "Illegal number of parameters"
        exit 1
    fi

    folder=$1
    analysis=$2
    num_archetypes=$3
    job_time=$4
    job_mem=$5
    job_cpus=$6
    
    # Check folder exists
    if [ ! -d "data/$folder" ] || [ ! -d "out/$folder" ]; then
        echo "Folder $folder does not exist"
        exit 1
    fi

    # Analysis can be either Archetypes or Seurat
    if [ "$analysis" != "Archetypes" ] && [ "$analysis" != "Seurat" ]; then
        echo "Analysis must be either Archetypes or Seurat"
        exit 1
    fi

    # Check num archetypes is a number >0
    if ! [[ "$num_archetypes" =~ ^[0-9]+$ ]]; then
        echo "Num archetypes must be a number"
        exit 1
    fi

    # Check job time is in the format HH:MM:SS
    if ! [[ "$job_time" =~ ^[0-9]+:[0-9]+:[0-9]+$ ]]; then
        echo "Job time must be in the format HH:MM:SS"
        exit 1
    fi

    # Check job mem is in the format <number>[K|M|G|T]
    if ! [[ "$job_mem" =~ ^[0-9]+[K|M|G|T]$ ]]; then
        echo "Job mem must be in the format <number>[K|M|G|T]"
        exit 1
    fi

    # Check job cpus is a number >0
    if ! [[ "$job_cpus" =~ ^[0-9]+$ ]]; then
        echo "Job cpus must be a number"
        exit 1
    fi
}

build_slurm() {
    local folder=$1
    local analysis=$2
    local num_archetypes=$3
    local job_time=$4
    local job_mem=$5
    local job_cpus=$6
    parse_parameters $folder $analysis $num_archetypes $job_time $job_mem $job_cpus

    local job_name="arch_$(basename $folder)"
    local output_dir="./out/${folder}"
    local slurm_output="${output_dir}/slurm_output.txt"
    local slurm_error="${output_dir}/slurm_error.txt"
    local singularity_log="${output_dir}/${job_name}.sing.log"
    local script_path="jobs/factory/${job_name}.sh"

    # mkdir -p "$output_dir"

    # Generate the job script
    cat > "$script_path" <<- EOF
#!/bin/bash

#SBATCH --job-name=$job_name
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=$slurm_output
#SBATCH --error=$slurm_error
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$job_cpus
#SBATCH --time=$job_time
#SBATCH --mem=$job_mem

# Command
singularity exec \\
  --env DISABLE_AUTH=true \\
  --env PASSWORD=psw \\
  --bind .:/app \\
  containers/singularityArchetypes.sif \\
  /app/jobs/archetypes.sh \\
    $analysis \\
    $folder \\
    $num_archetypes \\
  &> $singularity_log
EOF

    # Submit the job
    # sbatch "$script_path"
}

# local folder=$1
# local analysis=$2
# local num_archetypes=$3
# local job_time=$4
# local job_mem=$5
# local job_cpus=$6
build_slurm "Melanoma" "Archetypes" 13 "5:00:00" "150G" 15
build_slurm "MyocardialInfarction" "Archetypes" 13 "32:00:00" "280G" 20  
build_slurm "MouseCortex" "Archetypes" 13 "5:00:00" "150G" 15 

build_slurm "AllonKleinLab/Experiment1" "Archetypes" 13 "20:00:00" "220G"  15
build_slurm "AllonKleinLab/Experiment2" "Archetypes" 13 "32:00:00" "280G" 20 
build_slurm "AllonKleinLab/Experiment3" "Archetypes" 13 "10:00:00" "150G" 20 