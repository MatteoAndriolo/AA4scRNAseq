#!/bin/bash

#SBATCH --job-name=rExp1
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/AllonKleinLab/Experiment1/_slurm/%j.out.txt
#SBATCH --error=./out/AllonKleinLab/Experiment1/_slurm/%j.err.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250G

container="containers/img/myrubuntu.sif"
command="/app/Rmd/_main.sh -c Exp1 -w 20 --mink 11 --maxk 15 $@"

singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    $command \
  #&> ./out/AllonKleinLab/Experiment1/singularity.log
