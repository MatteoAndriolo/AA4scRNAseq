#!/bin/bash

#SBATCH --job-name=arch_Experiment3
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/AllonKleinLab/Experiment3/slurm_output.txt
#SBATCH --error=./out/AllonKleinLab/Experiment3/slurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=150G

# Command
singularity exec \
  --bind .:/app \
  containers/img/ubuntuArchetypes.sif \
  /app/singularity/AA.sh \
    Archetypes \
    AllonKleinLab/Experiment3 \
    10 \
  &> ./out/AllonKleinLab/Experiment3/singularity.log
#  --env DISABLE_AUTH=true \
#  --env PASSWORD=psw \
