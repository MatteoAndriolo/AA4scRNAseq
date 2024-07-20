#!/bin/bash

#SBATCH --job-name=rMelanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/_slurm/%j.out.txt
#SBATCH --error=./out/Melanoma/_slurm/%j.err.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=50G

container="containers/img/myrubuntu.sif"
command="/app/Rmd/_main.sh -c Melanoma -w 20 --mink 6 --maxk 18 $@ "

singularity exec \
  --bind .:/app \
    $container \
    $command 
