#!/bin/bash

#SBATCH --job-name=arch_Experiment2
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/AllonKleinLab/Experiment2/slurm_output.txt
#SBATCH --error=./out/AllonKleinLab/Experiment2/slurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=32:00:00
#SBATCH --mem=280G

# Command
singularity exec \
  --bind .:/app \
  containers/img/ubuntuArchetypes.sif \
  /app/singularity/AA.sh \
    Archetypes \
    AllonKleinLab/Experiment2 \
    9 \
  &> ./out/AllonKleinLab/Experiment2/singularity.log
#  --env DISABLE_AUTH=true \
#  --env PASSWORD=psw \
