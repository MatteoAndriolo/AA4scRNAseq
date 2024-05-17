#!/bin/bash

#SBATCH --job-name=render_Melanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/RMDslurm_output.txt
#SBATCH --error=./out/Melanoma/RMDslurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=1:00:00
#SBATCH --mem=50G

#
container=containers/img/ubuntuArchetypes.sif 
markdown=$1
# Rscript -e "rmarkdown::render()"
echo $markdown
singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    Rscript -e "rmarkdown::render($markdown)" \
  &> ./out/Melanoma/RMDsingularity.log
