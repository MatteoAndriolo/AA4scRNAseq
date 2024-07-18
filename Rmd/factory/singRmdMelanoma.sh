#!/bin/bash

#SBATCH --job-name=rMelanoma
#SBATCH --mail-user=nomatteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/_slurm/%j.out.txt
#SBATCH --error=./out/Melanoma/_slurm/%j.err.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=3:00:00
#SBATCH --mem=50G

container="containers/img/myrubuntu.sif"
command="/app/Rmd/_main.sh -c Melanoma $@ "
command="/app/Rmd/_main.sh -c Melanoma --mink 3 --maxk 10 $@ "

singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    $command \
  #&> ./out/Melanoma/${SLURM_JOB_ID}.log
