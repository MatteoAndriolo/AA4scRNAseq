#!/bin/bash

#SBATCH --job-name=AAMouse
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/MouseCortex/out.%t.txt
#SBATCH --error=./out/MouseCortex/err.%t.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G               
#SBATCH --time=15:00:00

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
 singularity exec  --bind .:/app $container  ${commands[i]} > ./out/MouseCortex/Par/out.$i.txt  &
 sleep 0.3
done

# Wait for all tasks to complete
wait
