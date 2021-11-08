#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Subsitute missing by -1
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
h5n <- args[1]
setPrefix <- gsub("assays.h5", "", h5n)

## Load libraries
library(DelayedMatrixStats)
library(HDF5Array)
library(SummarizedExperiment)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

assay(SE)[is.na(assay(SE))] <- -1

saveHDF5SummarizedExperiment(SE, ".", prefix = paste0(setPrefix, "missingSub_"))
