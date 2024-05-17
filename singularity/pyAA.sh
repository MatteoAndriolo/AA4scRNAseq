#!/bin/bash

#SBATCH --job-name=pyarch_Melanoma
#SBATCH --mail-user=matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type=ALL
#SBATCH --output=./out/Melanoma/slurm_output.txt
#SBATCH --error=./out/Melanoma/slurm_error.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=5:00:00
#SBATCH --mem=150G

# Command
singularity exec \
  --env DISABLE_AUTH=true \
  --env PASSWORD=psw \
  --bind .:/app \
  containers/singularityArchetypes.sif \
  /app/jobs/archetypes.sh \
    Archetypes \
    AllonKleinLab/Experiment1 \
    13 \
  &> ./out/AllonKleinLab/Experiment1/arch_Experiment1.sing.log
