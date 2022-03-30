#'#################################################################################
#'#################################################################################
#' Explore GSE57945 features from different models
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(HDF5Array)
library(pheatmap)
library(rjson)



makeHeatmap <- function(seObj){
  col_colors <- list(
    sample_type = c("Primary Tumor" = "lightgreen", "Solid Tissue Normal" = "black")
    # sex = c("Female" = "purple", "Male" = "lightblue")
  )

  pheatmap(assay(seObj), scale = "row",
           annotation_col  = data.frame(colData(seObj)[, c("sample_type"), drop = FALSE]),
           annotation_colors =  col_colors,
          show_colnames = FALSE)

}

kegg.map <- read.table("results/preprocess/kegg_filt_manual_gene_map.tsv", header = TRUE)
paths <- read.table("results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.2/model_trained/pathways_names.txt", header = TRUE)
genes <- read.table("results/TCGA_gexp_combat_coding/input_genes.txt", header = FALSE)

vst <- loadHDF5SummarizedExperiment("results/TCGA_gexp_combat_coding/", prefix = "vsd_norm")

path0 <- h5read("results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/relevance/relevance_genes_pathways.h5","patwhays_0")
rownames(path0) <- genes$V1
colnames(path0) <- colnames(vst)

genes0 <- subset(kegg.map, PathwayID == paths$X0[1])$Symbol

path0_sig <- path0[rownames(path0) %in% genes0, ]
pc0 <- prcomp(path0_sig)
