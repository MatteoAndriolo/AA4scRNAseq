#!/bin/bash

#SBATCH --job-name=rMelanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/_slurm/F%j.out.txt
#SBATCH --error=./out/Melanoma/_slurm/F%j.err.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=30:00:00
#SBATCH --mem=100G

container="containers/img/myrubuntu.sif"
command="/app/Rmd/_precomp.sh"

echo $container
echo $command

singularity exec \
  --bind .:/app \
    $container \
    $command 
