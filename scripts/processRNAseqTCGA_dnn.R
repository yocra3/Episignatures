#'#################################################################################
#'#################################################################################
#' Process TCGA gene expression for deep learning
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(topGO)
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)
register(MulticoreParam(10))

## Load data
load("data/tcga_gexp.Rdata")
gexp_fold <- "results/TCGA_gexp/"

dds <- DESeqDataSetFromMatrix(countData = assay(gexp_tcga),
                              colData = colData(gexp_tcga),
                              design = ~ project_id)
# dds <- DESeq(dds)
## Select genes expressed in at least 20 samples
keep <- rowSums(counts(dds) > 0) >= 20
vsd <- vst(dds[keep, ], blind=FALSE)

## Output vsd values
saveHDF5SummarizedExperiment(vsd, "results/TCGA_gexp/", prefix = "vsd_norm")


## Write gene names
genes <- rownames(vsd)
write.table(genes, file =  paste0(gexp_fold, "input_genes.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)


## Create labels file
project <- gexp_tcga$project_id
project[gexp_tcga$sample_type == "Solid Tissue Normal"] <- "Normal"


write.table(project, file = paste0(gexp_fold, "individuals_labels.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

## Select genes present in GOs
gexp_fold_go <- "results/TCGA_gexp_go/"

genevec <- factor(rep(c(0, 1), each = length(genes)/2))
names(genevec) <- genes

goData <- new("topGOdata",
              ontology = "BP",
              allGenes = as.factor(genevec),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "Ensembl")
vsd.go <- vsd[feasible(goData), ]

## Output vsd values
saveHDF5SummarizedExperiment(vsd.go, "results/TCGA_gexp_go/", prefix = "vsd_norm")


## Write gene names
genes <- rownames(vsd.go)
write.table(genes, file =  paste0(gexp_fold_go, "input_genes.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)


## Create labels file
project <- gexp_tcga$project_id
project[gexp_tcga$sample_type == "Solid Tissue Normal"] <- "Normal"
write.table(project, file = paste0(gexp_fold_go, "individuals_labels.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

## Train for association with sex
sex <- gexp_tcga$gender
sex[is.na(sex)] <- "female"
project2 <- paste(project, sex)
write.table(project2, file = paste0(gexp_fold_go, "individuals_labels_sex.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)


## Select autosomic
rowRanges(vsd.go) <- rowRanges(gexp_tcga)[rownames(vsd.go)]

vsd.go.auto <- vsd.go[!seqnames(vsd.go) %in% c("chrX", "chrY"), ]
saveHDF5SummarizedExperiment(vsd.go.auto, "results/TCGA_gexp_go/", prefix = "vsd_norm_auto")

genes.auto <- rownames(vsd.go.auto)
write.table(genes.auto, file =  paste0(gexp_fold_go, "input_genes_autosomics.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)
