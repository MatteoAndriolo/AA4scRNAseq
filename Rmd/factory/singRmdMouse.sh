#!/bin/bash

#SBATCH --job-name=rMouse
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/MouseCortex/_slurm/%j.out.txt
#SBATCH --error=./out/MouseCortex/_slurm/%j.err.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=150G

container="containers/img/myrubuntu.sif"
command="/app/Rmd/_main.sh -c Mouse -w 20 --mink 6 --maxk 18 $@"

singularity exec \
  --bind .:/app \
    $container \
    $command 
