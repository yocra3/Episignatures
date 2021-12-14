#'#################################################################################
#'#################################################################################
#' Correlation expression layer1 values
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(edgeR)
library(HDF5Array)

## Load layer values
tcga3_fold <- "results/TCGA_bioDNN/2021-11-21/model_features/bioDNN_v3/"
tcga3_conc <- read_delim(paste0(tcga3_fold, "concatenate.tsv"), delim = "\t")
tcga3_conc <- data.matrix(tcga3_conc)

labels <- read_delim("results/DNN/2021-11-19/model_trained/v3/layer_names.tsv", delim = "\t")
labels <- unlist(labels)

genelabels <- labels[!grepl("input|chr|dense", labels)]
genelabels <- genelabels[genelabels != "concatenate"]
colnames(tcga3_conc) <- c(genelabels, rep(paste0("chr", 1:22), each = 1000))

meth_se <- loadHDF5SummarizedExperiment("data/tcga_hdf5", prefix = "methyTCGA_")
rownames(tcga3_conc) <- colnames(meth_se)
rownames(tcga3_conc) <- substring(rownames(tcga3_conc), 1, 16)


tcga3_filt <- tcga3_conc[, !apply(tcga3_conc, 2, function(x) all(x == 0))]

## Load expression ####
load("data/tcga_gexp.Rdata")

dge <- DGEList(assay(gexp_tcga))
cpms <- cpm(dge)
colnames(cpms) <- substring(colnames(cpms), 1, 16)

# d0 <- calcNormFactors(d0)
# mm <- model.matrix(~1, colData(gexp_tcga))
# y <- voom(d0, mm, plot = T)

## Get common samples and "genes" ####
com_samps <- intersect(colnames(cpms), rownames(tcga3_filt))

tcga3_filt2 <- tcga3_filt[com_samps, ]
tcga3_filt2 <- tcga3_filt2[, !apply(tcga3_filt2, 2, function(x) all(x == 0))]

gexp_genes <- rowData(gexp_tcga)$external_gene_name
com_genes <- intersect(gexp_genes, colnames(tcga3_filt2))

### Remove duplicated genes
filt_genes <- gexp_genes[gexp_genes %in% com_genes]
dup_genes <- filt_genes[duplicated(filt_genes)]

com_genes <- com_genes[!com_genes %in% dup_genes]

cpms_filt <- cpms[gexp_genes %in% com_genes, com_samps]
rownames(cpms_filt) <- gexp_genes[gexp_genes %in% com_genes]

colnames(gexp_tcga) <- substring(colnames(gexp_tcga), 1, 16)
gexp_filt <- gexp_tcga[, com_samps]

meth_layer <- tcga3_filt2[, com_genes]

## Compute correlations ####
cors <- sapply(com_genes, function(g) cor(meth_layer[, g], cpms_filt[g, ]))

top_gene_cors <- names(cors)[!is.na(cors) & abs(cors) > 0.2]
all_cors <- cor(meth_layer, t(cpms_filt))
all_cors <- cor(meth_layer[, top_gene_cors], t(cpms_filt[top_gene_cors, ]))

cot <- apply(abs(all_cors), 1, rank)

pvals <- 1 - (diag(cot) -1)/ncol(cot)     

genes_tab <- data.frame(gene = colnames(cot), p_val = pvals, 
                        p_adj = p.adjust(pvals, method = "BH"),
                        r = diag(all_cors)) %>%
  mutate(sig = p_adj < 0.05) %>%
  as_tibble()
ggplot(genes_tab, aes(x = r, fill = sig)) + 
  geom_histogram() +
  facet_grid(~ sig)

plot(meth_layer[, "F12"], cpms_filt["F12", ])

## Compute prop activation per project
projects <- unique(gexp_filt$project_id)

meth_activation <- sapply(projects, function(proj) {
    
    sel_samps <- colnames(gexp_filt)[gexp_filt$project_id == proj]
    colMeans(meth_layer[sel_samps, ] > 0)
  }) 
  
gexp_activation <- sapply(projects, function(proj) {
  
  sel_samps <- colnames(gexp_filt)[gexp_filt$project_id == proj]
  rowMeans(cpms_filt[, sel_samps] > 1)
})
cors_act <- sapply(com_genes, function(g) cor(meth_activation[g, ], gexp_activation[g, ]))
plot(meth_layer[, "BRCA1"], cpms_filt["BRCA1", ])

corslist <- lapply(com_genes, function(g) {
  sapply(projects, function(proj) {
    
    sel_samps <- colnames(gexp_filt)[gexp_filt$project_id == proj]
    cor(meth_layer[sel_samps, g], cpms_filt[g, sel_samps])
  })
})


all_cor <- cor(t(cpms_filt), meth_layer)
