#!/bin/bash
dockername=$1   
name=$(basename $dockername)
noext=${name%.*}
echo singularity build $noext.sif docker-archive:$dockername
