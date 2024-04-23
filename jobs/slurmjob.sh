#!/bin/bash

#SBATCH --job-name archetypes
#SBATCH --mail-user matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type ALL

#SBATCH --output ./out/slurm_output.txt
#SBATCH --error ./out/slurm_error.txt

#SBATCH --partition allgroups
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --time 00:40
#SBATCH --mem 20G

# Necessary for Script
COMMAND="singularity exec /app/jobs/archetypes.sh"
# Execute command
sbatch $COMMAND