#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Select autosomic probes
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
setPrefix <- args[1]

## Load libraries
library(DelayedMatrixStats)
library(HDF5Array)
library(SummarizedExperiment)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

SE <- SE[seqnames(rowRanges(SE)) %in% c(1:22, paste0("chr", 1:22)), ]

saveHDF5SummarizedExperiment(SE, ".", prefix = paste0(setPrefix, "autosomicProbes_"))
