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

nworkers <- parallel::detectCores() - 2
plan("multicore", workers = nworkers)
