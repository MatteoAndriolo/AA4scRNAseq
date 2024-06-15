# Load necessary libraries
library(Seurat)
library(ggplot2)
library(here)
library(archetypes)
library(tidyr)
library(dplyr)
library(kneedle)
library(cowplot)
library(future)
# fetch number cpu
nworkers <- parallel::detectCores()
plan("multicore", workers = nworkers)
