# Load necessary libraries
library(Seurat)
library(archetypal)
library(archetypes)
library(cowplot)
library(parallel)
library(doParallel)
library(dplyr)
library(foreach)
library(future)
library(ggplot2)
library(here)
library(kneedle)
library(tidyr)

pathways <- list(
  "Glycolysis / Gluconeogenesis",
  "MAPK signaling pathway"
)

pathways <- list(
  "Glycolysis / Gluconeogenesis",
  "MAPK signaling pathway",
  "mTOR signaling pathway",
  "Pathways in cancer",
  "TGF-beta signaling pathway"
)
