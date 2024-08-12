#!/bin/bash

#SBATCH --job-name=AAMouse
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/MouseCortex/_slurm/%j.%t.out.txt
#SBATCH --error=./out/MouseCortex/_slurm/%j.%t.err.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=5              
#SBATCH --cpus-per-task=4       
#SBATCH --mem=200G               
#SBATCH --time=30:00:00

container="containers/img/myrubuntu.sif"

# Array of commands with different parameters
command="/app/Rmd/_main.sh -c Mouse -w 4 --mink 6 --maxk 18"
commands=(
  "$command -H TRUE"
  "$command -p 1"
  "$command -p 2"
  "$command -p 3"
  "$command -p 4"
  "$command -p 5"
)

# Loop through the commands array and run each one as a separate task using srun
for i in "${!commands[@]}"; do
  srun --exclusive --ntasks=1 \
    singularity exec \
        --bind ./app \
        $container \
        ${commands[i]} \
    &
done

# Wait for all tasks to complete
wait
