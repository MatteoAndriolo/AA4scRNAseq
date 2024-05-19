#!/bin/bash
command="singularity/singRmd.sh"
container="containers/img/ubuntuArchetypesPandas.sif"
sbatch $command $container md/Exp1.Rmd 
sbatch $command $container md/Exp2.Rmd 
sbatch $command $container md/Exp3.Rmd
sbatch $command $container md/Melanoma.Rmd