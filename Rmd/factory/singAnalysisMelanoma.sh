#!/bin/bash

#SBATCH --job-name=rMelanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/RMDslurm_output.txt
#SBATCH --error=./out/Melanoma/RMDslurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=1:00:00
#SBATCH --mem=100G

container="containers/img/myrubuntu.sif"


singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
    $container \
    $command \
  &> ./out/Melanoma/RMDsingularity.log

###
FOLDER_PATH="out/Melanoma/Melanoma_files/bigMelanomaHunique"
INPUT_FILE="MelanomaH_unique.rds"
OUTPUT_FILE="MelanomaH_analysis.rds"
command="_main_analysis.sh $FOLDER_PATH $INPUT_FILE $OUTPUT_FILE"

# Call the run_analysis.sh script with command line arguments
singularity exec \
  --bind .:/app \
  $container \
  $command \
  &> ./out/Melanoma/RMDanalysis.log