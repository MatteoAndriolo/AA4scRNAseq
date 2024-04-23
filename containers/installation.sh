#!/bin/bash

apt install -y libhdf5-dev libglpk-dev;

# Install R packages
#RUN R -e 'usethis::gitcreds_set(url = "https://github.com")'

R -e "\
    update.packages(ask = FALSE, checkBuilt = TRUE);\
    install.packages(c('usethis','rhdf5','remotes','devtools','BiocManager'));\
    install.packages('archetypes');\ 
    install.packages(c('units','sf','spdep'));\
    BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', \
            'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', \
            'batchelor', 'HDF5Array', 'terra', 'ggrastr'));\
    remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE);\
    remotes::install_github('satijalab/seurat-data', 'seurat5', quiet = TRUE);\
    remotes::install_github('satijalab/seurat-wrappers', 'seurat5', quiet = TRUE);\
    remotes::install_github('bnprks/BPCells', quiet = TRUE);\
    devtools::install_github('cole-trapnell-lab/monocle3');\
"


# RUN R -e 'remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE);'
# RUN R -e 'remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE);'
# RUN R -e 'setRepositories(ind=1:3) ; install.packages("Signac") '