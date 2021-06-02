#'#################################################################################
#'#################################################################################
#' Configuration for GSE109381 QC of methylation data
#'#################################################################################

## Discard samples for the project
discardSamples <- function(samplesheet){
  samplesheet 
}

## Create common sample ID between IDATs and phenotypes
addSampID <- function(samplesheet) {
  samplesheet$SampleID <- as.character(as.numeric(substring(samplesheet$Sample_Name, 7, 10)))
  samplesheet
}


## Variables from phenotype to check batch
batch_var <- c("Slide", "Array",  "Scan_Date", "Sample_Well", "Sample_Plate", "Sex", 
               "sges", "preterm", "BW", "BL", "HC", "tippart",
               "breastfeeding",  "msmk", "meduc", "edadm", "m_not_eur")

## QC parameters 
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1, 
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
  snp.concordance.threshold             = 0.95,
  sample.genotype.concordance.threshold = 0.8
)
pcs <- 20

## Annotation
array <- "IlluminaHumanMethylation450k"
annotation <- "ilmn12.hg19"