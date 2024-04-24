#!/bin/bash

#SBATCH --job-name arch_me
#SBATCH --mail-user matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type ALL

#SBATCH --output ./out/Melanoma/slurm_output.txt
#SBATCH --error ./out/Melanoma/slurm_error.txt

#SBATCH --partition allgroups
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --time 40:00
#SBATCH --mem 30G

ANALYSIS=Archetypes
DATASET=Melanoma
NUM_ARCHETYPES=13
SLURM_LOG=out/$DATASET/slurm.log

WORKSPACE="/home/andrioloma/MasterThesis"

# COMMAND
singularity exec \
    --env DISABLE_AUTH=true \
    --env PASSWORD=psw \
    --bind $WORKSPACE:/app \
     $WORKSPACE/containers/singularityArchetypes.sif \
     /app/jobs/archetypes.sh \
        $ANALYSIS \
        $DATASET \
        $NUM_ARCHETYPES \
     &> $WORKSPACE/out/$DATASET/$ANALYSIS.$NUM_ARCHETYPES.sing.log


