#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Remove probes with methylation range < 0.1
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

quant <- rowQuantiles(assay(SE), probs = c(0.99, 0.01), na.rm = TRUE)
ranges <- apply(quant, 1, function(x) x[1] - x[2])
SE <- SE[ranges > 0.1, ]

saveHDF5SummarizedExperiment(SE, ".", prefix = paste0(setPrefix, "variantProbes_"))
