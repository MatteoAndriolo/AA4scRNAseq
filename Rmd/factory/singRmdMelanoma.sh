#!/bin/bash

#SBATCH --job-name=rMelanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/RMDslurm_output.txt
#SBATCH --error=./out/Melanoma/RMDslurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2:00:00
#SBATCH --mem=100G

container="containers/img/ubuntuArchetypesPandas.sif"
container="containers/img/myrocker.sif"
container="containers/img/myrubuntu.sif"

# command="Rscript -e \"rmarkdown::render('Rmd/Melanoma.Rmd')\""
command="/app/Rmd/_main.sh -c Melanoma $@"
#command="/app/Rmd/_main.sh -c Melanoma $@"

singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    $command \
  &> ./out/Melanoma/RMDsingularity.log
