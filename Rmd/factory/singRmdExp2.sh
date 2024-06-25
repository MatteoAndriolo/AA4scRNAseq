#!/bin/bash

#SBATCH --job-name=rExp2
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/AllonKleinLab/Experiment2/RMDslurm_output.txt
#SBATCH --error=./out/AllonKleinLab/Experiment2/RMDslurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=4:00:00
#SBATCH --mem=200G

container="containers/img/ubuntuArchetypesPandas.sif"
container="containers/img/myrocker.sif"
container="containers/img/myrubuntu.sif"

#command="Rscript -e \"rmarkdown::render('Rmd/Exp2.Rmd')\"" 
command="/app/Rmd/_main.sh -c Exp2 $@"

singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    $command \
  &> ./out/AllonKleinLab/Experiment2/RMDsingularity.log
