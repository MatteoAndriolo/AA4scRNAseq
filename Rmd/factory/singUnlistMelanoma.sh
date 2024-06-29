#!/bin/bash

#SBATCH --job-name=rAnMelanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/RMDUnSlurm_output.txt
#SBATCH --error=./out/Melanoma/RMDUnSlurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=20:00
#SBATCH --mem=50G


FOLDER_PATH="out/Melanoma/Melanoma_files/bigMelanomaHunique"
INPUT_FILE="MelanomaH_unique.rds"
OUTPUT_FILE="MelanomaH_analysis.rds"

command="Rscript /app/Rmd/unclassObject.R"
container="containers/img/myrubuntu.sif"

# Call the run_analysis.sh script with command line arguments
singularity exec \
  --bind .:/app \
  $container \
  $command 