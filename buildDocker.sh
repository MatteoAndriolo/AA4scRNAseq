#!/bin/sh

if [ "$1" = "rocker" ]; then
    docker buildx build -t myrocker:latest --build-arg GITHUB_PAT=$(cat containers/githubpat) -f containers/Dockerfile/Dockerfile.AA.rocker .
    docker save -o cImages/myrocker.tar myrocker:latest
    singularity build myrocker.sif docker-archive://cImages/myrocker.tar
elif [ "$1" = "ubuntu" ]; then
    docker buildx build -t myrubuntu:latest --build-arg GITHUB_PAT=$(cat containers/githubpat) -f containers/Dockerfile/Dockerfile.AA.ubuntu .
    docker save -o cImages/myrubuntu.tar myrubuntu:latest
    singularity build myrubuntu.sif docker-archive://cImages/myrubuntu.tar
else 
    echo "Usage: ./buildDocker.sh ubuntu|rocker"
fi
