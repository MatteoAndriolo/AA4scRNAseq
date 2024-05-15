#!/bin/bash

#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --partition=allgroups
#SBATCH --job-name=test_singularity
#SBATCH --output=result-%j.out
#SBATCH --error=error-%j.err
#SBATCH --time=00:01:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

container="containers/img/ubuntuArchetypes.sif"

# Your Singularity command here
singularity exec --bind .:/app $container env > enctext

