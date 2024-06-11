#!/bin/bash

#SBATCH --job-name=rMelanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/RMDslurm_output.txt
#SBATCH --error=./out/Melanoma/RMDslurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=10:00:00
#SBATCH --mem=150G

container="containers/img/ubuntuArchetypesPandas.sif"

# command="Rscript -e \"rmarkdown::render('Rmd/Melanoma.Rmd')\""
command="/app/Rmd/_main.sh Melanoma"

singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    $command \
  &> ./out/Melanoma/RMDsingularity.log
