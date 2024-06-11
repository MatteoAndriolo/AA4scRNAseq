# Load necessary libraries
library(Seurat)
library(ggplot2)
library(here)
library(archetypes)
library(tidyr)
library(dplyr)
if (!require(kneedle)) {
  install.packages("quantmod")
  devtools::install_github("etam4260/kneedle")
}
library(kneedle)
library(cowplot)
