#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Output projects labels
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
setPrefix <- args[1]

## Load libraries
library(HDF5Array)
library(SummarizedExperiment)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

project <- SE$project
project[SE$sample_type == "Solid Tissue Normal"] <- "Normal"

write.table(project, file = "TCGA_individuals_cancer_labels.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

