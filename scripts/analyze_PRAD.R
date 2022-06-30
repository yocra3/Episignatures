#'#################################################################################
#'#################################################################################
#' Analyze PRAD
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(fgsea)
library(e1071)
library(HDF5Array)
library(hipathia)
library(org.Hs.eg.db)



load("data/tcga_gexp_combat.Rdata")
genes <- read.table("./results/GTEx_coding/input_genes.txt")

path.map <- read.table("results/GTEx_coding/go_kegg_final_gene_map.tsv", header = TRUE)

prad.feat <- read.table("results/GTEx_coding_PRAD_tumor/paths_filt2_full_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(prad.feat) <- paths.vec

## Subset data
prad.vst <- loadHDF5SummarizedExperiment("results/TCGA_gexp_coding_noPRAD/", prefix = "vsd_norm_prad_tumor")

prad <- gexp_tcga_combat[genes$V1, gexp_tcga_combat$project_id == "TCGA-PRAD"]
prad <- prad[, !is.na(prad$paper_Reviewed_Gleason_category)]

## DE in raw data
prad$gleason <- ifelse(prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low")
ddsSE <- DESeqDataSet(prad, design = ~ paper_Subtype + age_at_index + race + gleason )
vst.prad <- vst(ddsSE, blind=FALSE)

pc.ori <- prcomp(t(assay(vst.prad)))
data.frame(pc.ori$x, var = prad$paper_Reviewed_Gleason_category) %>% ggplot(aes(x = PC1, y = PC2, color = var)) + geom_point()
data.frame(pc.ori$x, var = prad$gleason) %>% ggplot(aes(x = PC1, y = PC2, color = var)) + geom_point()
data.frame(pc.ori$x, var = prad$paper_Subtype) %>% ggplot(aes(x = PC1, y = PC2, color = var)) + geom_point()


dds <- DESeq(ddsSE)
res_prad <- results(dds)
res_prad$p.adj.bf <- p.adjust(res_prad$pvalue )
save(res_prad, file = "results/TCGA_PRAD/genes_results.Rdata")

## Run GO enrichment
# dif.genes <- as.numeric(res$padj < 0.05)
# dif.genes[is.na(dif.genes)] <- 0
# names(dif.genes) <- rownames(res)
#
# pwf <- nullp(dif.genes, "hg19", "ensGene")
# gos <- goseq(pwf, "hg19", "ensGene", test.cats = c("GO:BP"))
# gos.mod <- subset(gos, category %in% paths.vec)
# gos.mod$prop_DE <- gos.mod$numDEInCat / gos.mod$numInCat
# rownames(gos.mod) <- gos.mod$category

## Run pathways differential Analysis
mod <- model.matrix(~ gleason + paper_Subtype + age_at_index + race, colData(prad))
lm.path <- lmFit(t(prad.feat), mod) %>% eBayes()
tab.path_prad <- topTable(lm.path, coef = 2, n = Inf)
tab.path_prad$category <- rownames(tab.path_prad)

tab.path_prad$DE_prop <- sapply( tab.path_prad$category, function(cat) {
  genes <- subset(path.map, PathwayID  == cat )$Symbol
  mini_tab <- res_prad[genes, ]
  mean(mini_tab$padj  < 0.05, na.rm = TRUE)
})
save(tab.path_prad, file = "results/TCGA_PRAD/pathways_results.Rdata")


png("figures/TCGA_propDE_vs_pvalPaths.png")
ggplot(tab.path_prad, aes(x = DE_prop, y = -log10(P.Value ))) +
 geom_point() +
 scale_x_continuous(name = "Proportion of genes DE") +
 scale_y_continuous(name = "-log10 P-value Pathways") +
 theme_bw()
dev.off()
cor(tab.path_prad$ DE_prop, -log10(tab.path_prad$P.Value), use = "complete")
# [1] 0.3421541

cor(tab.path_prad$DE_prop, abs(tab.path_prad$logFC), use = "complete")
# 0.327146


cor(tab.path_prad$DE_prop[tab.path_prad$DE_prop > 0 ], abs(tab.path_prad$logFC)[tab.path_prad$DE_prop > 0 ], use = "complete")
# 0.3294506
# comb.df <- left_join(tab.path, gos.mod, by = "category")
# comb.df$min_p <- pmin(comb.df$over_represented_pvalue, comb.df$under_represented_pvalue)
# plot(-log10(comb.df$min_p), -log10(comb.df$P.Value ))
# cor(-log10(comb.df$min_p), -log10(comb.df$P.Value ), use = "complete")

## fgsea
ranks <- res_prad$log2FoldChange
ranks[is.na(ranks)] <- 0

entrez_ids <- mapIds(org.Hs.eg.db, rownames(res_prad), keytype="ENSEMBL", column="ENTREZID")
names(ranks) <- entrez_ids

pathways <- reactomePathways(entrez_ids)
# names(paths.vec) <- paths.vec
# pathways <- lapply( paths.vec, function(x) subset(path.map, PathwayID  == x )$Symbol)
fgseaRes_prad <- fgsea(pathways = pathways, stats = ranks)
save(fgseaRes_prad, file = "results/TCGA_PRAD/fgsea_results.Rdata")
#
# df.comb2 <- left_join(fgseaRes, tab.path, by = "category")
#
# png("figures/TCGA_pval_gseavsPaths.png")
# ggplot(df.comb2, aes(x = -log10(P.Value ), y = -log10(pval ))) +
#  geom_point() +
#  scale_y_continuous(name = "-log10 P-value GSEA") +
#  scale_x_continuous(name = "-log10 P-value Pathways") +
#  theme_bw()
# dev.off()


## hipathia
trans_data <- translate_data(vst.prad, "hsa")
exp_data <- normalize_data(trans_data)
hip_pathways <- load_pathways(species = "hsa")

hip.res_prad <- hipathia(exp_data, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.path_vals <- get_paths_data(hip.res_prad )
hip.comp_prad <- do_wilcoxon(hip.path_vals, hip.path_vals$gleason, g1 = "High", g2 = "Low")

save(hip.comp_prad, file = "results/TCGA_PRAD/hipathia.res.Rdata")


## Train SVM to classify gleason
colnames(prad.feat) <- gsub(":", "_", colnames( prad.feat))

### All features
df_svm <-  data.frame(gleason = factor(prad$gleason), prad.feat)
svm_gleason <- svm(gleason ~ ., df_svm)
pred.tcga <- predict(svm_gleason, prad.feat)
table(prediction = pred.tcga , real = prad$gleason )

sel_paths <- subset(tab.path_prad, adj.P.Val < 0.05)
svm_gleason_filt <- svm(gleason ~ ., df_svm[, c("gleason", gsub(":", "_", rownames( sel_paths)))])
pred.tcga_filt <- predict(svm_gleason_filt, prad.feat)
table(prediction = pred.tcga_filt , real = prad$gleason )



save(svm_gleason, svm_gleason_filt, file = "results/TCGA_PRAD/svm_model.Rdata")

load("results/TCGA_PRAD/svm_model_geo.Rdata")
pred.tcga2 <- predict(svm_gleason_geo, prad.feat)
table(prediction = pred.tcga2 , real = prad$gleason )

pred.tcga_filt2 <- predict(svm_gleason_geo_filt, prad.feat)
table(prediction = pred.tcga_filt2 , real = prad$gleason )
