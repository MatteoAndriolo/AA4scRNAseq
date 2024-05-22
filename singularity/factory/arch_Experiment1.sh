#!/bin/bash

#SBATCH --job-name=arch_Experiment1
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/AllonKleinLab/Experiment1/slurm_output.txt
#SBATCH --error=./out/AllonKleinLab/Experiment1/slurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=20:00:00
#SBATCH --mem=220G

# Command
singularity exec \
  --bind .:/app \
  containers/img/ubuntuArchetypes.sif \
  /app/singularity/AA.sh \
    Archetypes \
    AllonKleinLab/Experiment1 \
    8 \
  &> ./out/AllonKleinLab/Experiment1/singularity.log
#  --env DISABLE_AUTH=true \
#  --env PASSWORD=psw \
