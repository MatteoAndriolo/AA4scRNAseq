install.packages(c("remotes","devtools", "rhdf5"))
usethis::gitcreds_set()

library(remotes)
remotes::git_credentials("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")
update.packages(ask = FALSE, checkBuilt = TRUE)

install.packages("archetypes")

remotes::install_github("satijalab/seurat", "seurat5")
remotes::install_github("satijalab/seurat-data", "seurat5")
remotes::install_github("satijalab/azimuth", "seurat5")
remotes::install_github("satijalab/seurat-wrappers", "seurat5")
remotes::install_github("stuart-lab/signac", "seurat5")

remotes::install_github("satijalab/seurat", "seurat5", quiet=TRUE)
remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
