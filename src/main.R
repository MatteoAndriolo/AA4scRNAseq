#!/usr/bin/env Rscript
source("src/archetypes.R")


# Function to display help information
displayHelp <- function() {
  cat("Help Information:\n")
  cat("Valid Analysis Types:\n")
  cat(paste("  -", paste(validAnalysisTypes, collapse="\n  -")), "\n\n")
  cat("Valid Datasets:\n")
  cat(paste("  -", paste(validDatasets, collapse="\n  -")), "\n\n")
  cat("Usage:\n")
  cat("  Rscript file.R <analysisType> <dataset> [<numArchetypes>]\n")
  cat("Example:\n")
  cat("  Rscript file.R Seurat MouseCortex\n")
  cat("  Rscript file.R Archetypes AllenKlein-Exp1 5\n")
}

# Main function to decide which analysis to run
mainAnalysisFunction <-
  function(analysisType, dataset, numArchetypes = NULL) {
    validAnalysis = c("Seurat", "Archetypes")
    if (!analysisType %in% validAnalysis) {
      stop("Invalid analysis type. Choose either 'Seurat' or 'Archetypes'.")
    }
    
    if (!dataset %in% c(
      # "AllonKleinLab-Exp1",
      # "AllonKleinLab-Exp2",
      # "AllonKleinLab-Exp3",
      "Melanoma",
      # "MouseCortex",
      "MyocardialInfarction"
    )) {
      stop("Invalid dataset name.")
    }
    
    
    if (getwd() != "/app") {
      setwd("/app")
    }
    getwd()
    runArchetypeAnalysis(dataset, numArchetypes)
  }


# Read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript file.R <analysisType> <dataset> [<numArchetypes>]")
}

analysisType <- args[1]
dataset <- args[2]
numArchetypes <-
  ifelse(length(args) >= 3, as.numeric(args[3]), NULL)

# Call main function with the provided arguments
mainAnalysisFunction(analysisType, dataset, numArchetypes)

