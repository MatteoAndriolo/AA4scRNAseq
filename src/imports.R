# Load necessary libraries
if (!require(networkD3)) install.packages("networkD3")
library(Rtsne)
library(Seurat)
library(archetypal)
#library(archetypes)
library(cowplot)
library(doParallel)
library(dplyr)
library(foreach)
library(future)
library(ggplot2)
if (!require(ggsankey)) devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(here)
if (!require(kneedle)) devtools::install_github("etam4260/kneedle")
library(kneedle)
library(parallel)
library(plotly)
library(scales)
library(stringr)
if (!require(tidyverse)) install.packages("tidyverse")
library(tidyr)
library(viridis)

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

namePathw <- list(
  "GLYK",
  "MAPK",
  "MTOR",
  "CANCER",
  "TGF",
  "HVF"
)
