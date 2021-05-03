#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Remove probes with call rate < 90%
#' Remove non-CpG probes
#' Sort object by genomic coordinates
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

## Remove non-cg positions
SE <- SE[grep("cg", rownames(SE)), ]

## Filter CpGs with all missings
pNA <- rowMeans(is.na(assay(SE)))
SE <- SE[pNA < 0.9, ]

## sort SE by GenomicCoordinates
SE <- sort(SE)
saveHDF5SummarizedExperiment(SE, ".", prefix = paste0(setPrefix, "probesFilt_"))

