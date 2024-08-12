#!/bin/bash

#SBATCH --job-name=rExp3
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/AllonKleinLab/Experiment3/_slurm/%j.out.txt
#SBATCH --error=./out/AllonKleinLab/Experiment3/_slurm/%j.err.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250G

container="containers/img/myrubuntu.sif"
command="/app/Rmd/_main.sh -c Exp3 -w 20 --mink 12 --maxk 17 $@"

singularity exec \
  --bind .:/app \
    $container \
    $command 
