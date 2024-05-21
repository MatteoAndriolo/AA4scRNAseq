#!/bin/bash

#SBATCH --job-name=rMelanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/RMDslurm_output.txt
#SBATCH --error=./out/Melanoma/RMDslurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=5:00:00
#SBATCH --mem=150G
container="containers/img/ubuntuArchetypesPandas.sif"
singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    Rscript -e "rmarkdown::render('md/Melanoma.Rmd')" \
  &> ./out/Melanoma/RMDsingularity.log
