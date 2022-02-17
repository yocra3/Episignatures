#'#################################################################################
#'#################################################################################
# Explore genes grouped in neurons
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(topGO)
library(tidyverse)
library(HDF5Array)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(parallel)


## Load data ####
genes <- read.table("results/TCGA_gexp/input_genes.txt")
genes <- as.character(genes$V1)

gen_mat <- read.table("results/TCGA_gexp/v1.4.prun.1/model_trained/selFeatures_dense.tsv",
  sep = ",", header = TRUE)
gen_mat <- data.matrix(gen_mat[, -1])
rownames(gen_mat) <- genes


genes <- read.table("results/TCGA_gexp_go/input_genes_autosomics.txt")
genes <- as.character(genes$V1)

gen_mat <- read.table(gzfile("results/TCGA_gexp_norm/v1.4.prun.2/model_trained/dense_weights.txt.gz"))
gen_mat <- data.matrix(gen_mat)
rownames(gen_mat) <- genes
colnames(gen_mat) <- paste0("X", 0:(ncol(gen_mat) - 1))

## Select columns with at least 1 gene
gen_mat_f <- gen_mat[, colSums(gen_mat) > 0]

## Convert matrix layer to list for topGO
topList <- lapply(seq_len(ncol(gen_mat_f)), function(i) {
  vec <- gen_mat_f[, i]
  names(vec) <- genes
  vec
})


goData <- new("topGOdata",
            ontology = "BP",
            allGenes = as.factor(topList[[1]]),
            annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "Ensembl")

computeGOs <- function(genesv){

  Data <- updateGenes(goData, as.factor(genesv))
  mixed <- runTest(Data, algorithm = "classic", statistic = "fisher")


  finTab <- GenTable(Data,
                    res = mixed,
                    topNodes = length(score(mixed)))
  finTab$n_genes <- sum(genesv)
  finTab$prop_genes_neuron <- finTab$Significant/finTab$n_genes
  finTab$prop_genes_GO <- finTab$Significant/finTab$Annotated

  finTab
}
gos_list <- mclapply(topList, computeGOs, mc.cores = 10)

go_genes <- allGenes(goData)[feasible(goData)]
gos_sumdf <- data.frame(propGO =  sapply(gos_list, function(x) {
                          x.f <- subset(x, prop_genes_neuron == max(x$prop_genes_neuron))
                          max(x.f$prop_genes_GO)
                        }),
                        found =  sapply(gos_list, function(x) {
                                                  x.f <- subset(x, prop_genes_neuron == max(x$prop_genes_neuron))
                                                  max(x.f$Significant)
                                                }),
                        nGO =  sapply(gos_list, function(x) {
                          x.f <- subset(x, prop_genes_neuron == max(x$prop_genes_neuron))
                          x.f2 <-  subset(x.f, prop_genes_GO == max(x.f$prop_genes_GO))
                          max(x.f2$Annotated)
                        } ),
                        nNeuron =  sapply(topList, function(x) {
                                                  sel <- names(x)[x == 1]
                                                  sel <- sel[sel %in% go_genes]
                                                  length(sel)
                                                })) %>%
      as_tibble()
gos_sumdf <- gos_sumdf %>%
  mutate(propNeuron = found/nNeuron)


## Convert logical layer to gene list
gene_list <- lapply(seq_len(ncol(gen_mat_f)), function(i) genes[as.logical(gen_mat_f[, i])])
