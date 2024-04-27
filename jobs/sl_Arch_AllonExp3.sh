#!/bin/bash

#SBATCH --job-name arch_me
#SBATCH --mail-user matteo.andriolo.2@studenti.unipd.it
#SBATCH --mail-type ALL

#SBATCH --output ./out/AllonKleinLab/Experiment3/slurm_output.txt
#SBATCH --error ./out/AllonKleinLab/Experiment3/slurm_error.txt

#SBATCH --partition allgroups
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 5:00:00
#SBATCH --mem 150GB

ANALYSIS=Archetypes
DATASET=AllonKleinLab/Experiment3
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


