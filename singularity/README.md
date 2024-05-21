# Singularity folder structure:

This folder contains bash script that will be runned by slurm usin `sbatch` command

## Files:

- `mainAA.sh` :
- `_slurmFactory` : used to generate slurm job descriptions
  - files are generated in `factory` folder
- `_AA.sh` : used by `_slurmFactory`. Runs Rscript for given parameters:
  - $SCRIPT_PATH
  - $ANALYSIS
  - $DATASET
  - $NUM_ARCHETYPES
- `main.sh` script that starts all jobs present in `factory` folder
