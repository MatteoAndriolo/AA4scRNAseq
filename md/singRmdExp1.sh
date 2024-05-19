#!/bin/bash

#SBATCH --job-name=rExp1
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/AllonKleinLab/Experiment1/RMDslurm_output.txt
#SBATCH --error=./out/AllonKleinLab/Experiment1/RMDslurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=5:00:00
#SBATCH --mem=200G

singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    Rscript -e "rmarkdown::render('md/Exp1.Rmd')" \
  &> ./out/AllonKleinLab/Experiment1/RMDsingularity.log
