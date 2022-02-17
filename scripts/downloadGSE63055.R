#'#################################################################################
#'#################################################################################
#' Process GSE63055 gene expression for deep learning
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(recount)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)
register(MulticoreParam(10))

## Download the RangedSummarizedExperiment object at the gene level for study SRP049593
url <- download_study('SRP049593', outdir = 'data/SRP049593')

## Load the data
load(file.path('data/SRP049593', 'rse_gene.Rdata'))


url <- download_study('SRP042228', outdir = 'data/SRP042228')
load(file.path('data/ERP009437', 'rse_gene.Rdata'))
