#!/bin/bash

#SBATCH --job-name=reductions
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=/home/andrioloma/MasterThesis/createImgOut.slurm
#SBATCH --error=/home/andrioloma/MasterThesis/createImgOut.slurm
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --mem=250G

command="Rscript src/VisualizeData.R"
#container="containers/img/singularityArchetypes.sif"
container="containers/img/ubuntuArchetypes.sif"

# Command
#singularity exec --env DISABLE_AUTH=true --env PASSWORD=psw --fakeroot --bind /home/andrioloma/MasterThesis:/app $container $command &> /home/andrioloma/MasterThesis/createImgOut.sing
singularity exec --bind .:/app $container $command &> /home/andrioloma/MasterThesis/createImgOut.sing
