#'#################################################################################
#'#################################################################################
#' Process TCGA gene expression for deep learning
#'#################################################################################
#'#################################################################################


## Load libraries ####
# library(topGO)
library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)
library(sva)
library(rtracklayer)

register(MulticoreParam(10))

## Load data
load("data/tcga_gexp.Rdata")
gexp_fold <- "results/TCGA_gexp_combat/"

### Correct sequencing center with ComBat
phenos <- get_IDs(gexp_tcga)


gexp_tcga$sample_type2 <- ifelse(gexp_tcga$sample_type == "Solid Tissue Normal", "normal", "tumor")

adj.mod <- model.matrix(~ project_id, colData(gexp_tcga))[, -14]## Avoid confounding
adj_counts <- ComBat_seq(assay(gexp_tcga), batch = phenos$center , group = gexp_tcga$sample_type2, covar_mod = adj.mod, full_mod=TRUE)

gexp_tcga_combat <- gexp_tcga
assay(gexp_tcga_combat) <- adj_counts
save(gexp_tcga_combat, file = "data/tcga_gexp_combat.Rdata")

## Load gencode annotation
annot <- readGFF("/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz")
annot <- filter(annot, type == "gene") %>%
  as_tibble()

annot_cod <- subset(annot, gene_type == "protein_coding")   %>%
 mutate(ensembl = gsub("\\.[0-9]*_[0-9]*", "", gene_id , perl = TRUE))

adj_counts2 <- adj_counts
mode(adj_counts2) <- "integer"
adj_counts2[is.na(adj_counts2)] <- max(adj_counts2, na.rm = TRUE)
dds <- DESeqDataSetFromMatrix(countData = adj_counts2,
                              colData = colData(gexp_tcga),
                              design = ~ project_id)
# dds <- DESeq(dds)
## Select genes expressed in at least 20 samples
keep <- rowSums(counts(dds) > 0) >= 20
vst <- vst(dds[keep, ], blind=FALSE)

## Output vsd values
saveHDF5SummarizedExperiment(vst, "results/TCGA_gexp_combat/", prefix = "vsd_norm")

## Write gene names
genes <- rownames(vst)
write.table(genes, file =  paste0(gexp_fold, "input_genes.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)


## Create labels file
project <- gexp_tcga$project_id
project[gexp_tcga$sample_type == "Solid Tissue Normal"] <- "Normal"


write.table(project, file = paste0(gexp_fold, "individuals_labels.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

## Select protein coding genes
annot <- readGFF("/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz")
annot <- filter(annot, type == "gene") %>%
  as_tibble() %>%
     mutate(Symbol = gsub("_[0-9]*", "", gene_id , perl = TRUE),
          Symbol = gsub("\\.[0-9]*", "", Symbol , perl = TRUE))

cod_genes <- subset(annot, gene_type == "protein_coding")$Symbol
vst_coding <- vst[rownames(vst) %in% cod_genes, ]

gexp_fold_cod <- "results/TCGA_gexp_combat_coding/"

saveHDF5SummarizedExperiment(vst_coding, "results/TCGA_gexp_combat_coding/", prefix = "vsd_norm")

## Write gene names
genes <- rownames(vst_coding)
write.table(genes, file =  paste0(gexp_fold_cod, "input_genes.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)


## Create labels file
project <- as.character(vst_coding$project_id)
project[vst_coding$sample_type == "Solid Tissue Normal"] <- "Normal"
write.table(project, file = paste0(gexp_fold_cod, "individuals_labels.txt"), quote = FALSE,
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
