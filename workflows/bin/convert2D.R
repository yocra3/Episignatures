#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Add distance and CpGs distribution to convert to 2D
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
setPrefix <- args[1]

## Load libraries
library(DelayedMatrixStats)
library(HDF5Array)
library(SummarizedExperiment)
library(rhdf5)

impute.matrix <- function (x, margin = 1, fun = function(x) mean(x, na.rm = T)) 
{
  if (margin == 2) 
    x <- t(x)
  idx <- which(is.na(x) | !is.finite(x), arr.ind = T)
  if (length(idx) > 0) {
    na.idx <- unique(idx[, "row"])
    v <- apply(x[na.idx, , drop = F], 1, fun)
    v[which(is.na(v))] <- fun(v)
    x[idx] <- v[match(idx[, "row"], na.idx)]
    stopifnot(all(!is.na(x)))
  }
  if (margin == 2) 
    x <- t(x)
  x
}
SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

## Prepare quantile
quant <- rowQuantiles(assay(SE), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)

assay(SE) <- impute.matrix(data.matrix(assay(SE)), margin = 2)
h5createFile("assay_conv_2D.h5")
h5createDataset("assay_conv_2D.h5", "methy", 
                dims = list(nrow(SE), ncol(quant) + 2, ncol(SE)))

for (x in seq_len(ncol(SE))){
  val <- cbind(assay(SE[, x]), quant, start(SE))
  h5write(val, "assay_conv_2D.h5", "methy", index = list(NULL, NULL, x))
} 
h5closeAll()



