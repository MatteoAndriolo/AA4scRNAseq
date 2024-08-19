#!/bin/bash

#SBATCH --job-name=AAplots
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=30:00
#SBATCH --mem=50G

container="containers/img/myrubuntu.sif"
command="Rscript /app/Rmd/z_analysisIterColors.R "

singularity exec \
  --bind .:/app \
    $container \
    $command 
