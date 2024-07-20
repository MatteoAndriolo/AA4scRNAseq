#!/bin/bash

#SBATCH --job-name=rExp2
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/AllonKleinLab/Experiment2/_slurm/%j.out.txt
#SBATCH --error=./out/AllonKleinLab/Experiment2/_slurm/%j.err.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250G

container="containers/img/myrubuntu.sif"
command="/app/Rmd/_main.sh -c Exp2 -w 20 --mink 10 --maxk 15 $@"

singularity exec \
  --bind .:/app \
    $container \
    $command 
