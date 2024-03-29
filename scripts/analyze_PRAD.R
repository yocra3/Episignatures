docker run -it -v /home/SHARED/PROJECTS/Episignatures:/home/SHARED/PROJECTS/Episignatures -w "$PWD" yocra3/episignatures_rsession:1.3  /bin/bash # nolint: error.
R
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
library(GSVA)
library(e1071)
library(HDF5Array)
library(hipathia)
library(org.Hs.eg.db)
library(parallel)
library(NetActivity)
library(NetActivityData)

load("data/tcga_gexp_combat.Rdata")

## Prepare PRAD data
prad.all <- gexp_tcga_combat[, gexp_tcga_combat$project_id == "TCGA-PRAD"]
prad <- prad.all[, !is.na(prad.all$paper_Reviewed_Gleason_category)]
prad$gleason <- factor(ifelse(prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low"), levels = c("Low", "High"))

ddsSE <- DESeqDataSet(prad, design = ~ paper_Subtype + age_at_index + race + gleason )
vst.prad <- vst(ddsSE, blind=FALSE)
save(vst.prad, file = "results/TCGA_PRAD/vst_SE.Rdata")


## DE in raw data
dds <- DESeq(ddsSE)
res_prad <- results(dds)
res_prad$p.adj.bf <- p.adjust(res_prad$pvalue )
save(res_prad, file = "results/TCGA_PRAD/genes_results.Rdata")

## NetActivity
preproc_prad <- prepareSummarizedExperiment(vst.prad, "gtex_gokegg")
scores_prad <- computeGeneSetScores(preproc_prad, "gtex_gokegg")
save(scores_prad, file = "results/TCGA_PRAD/NetActivity_scores.Rdata")


mod <- model.matrix(~ gleason + paper_Subtype + age_at_index + race, colData(scores_prad))
lm.path <- lmFit(assay(scores_prad), mod) %>% eBayes()
tab.path_prad <- topTable(lm.path, coef = 2, n = Inf)
tab.path_prad$category <- rownames(tab.path_prad)

tab.path_prad$DE_prop <- sapply( tab.path_prad$category, function(cat) {
  genes <- subset(path.map, PathwayID  == cat )$Symbol
  mini_tab <- res_prad[genes, ]
  mean(mini_tab$padj  < 0.05, na.rm = TRUE)
})

# tab.path_prad_race <- topTable(lm.path, coef = 13, n = Inf)
# tab.path_prad_race$category <- rownames(tab.path_prad_race)

save(tab.path_prad, file = "results/TCGA_PRAD/pathways_results.Rdata")


png("figures/TCGA_propDE_vs_pvalPaths.png")
ggplot(tab.path_prad, aes(x = DE_prop, y = -log10(P.Value ))) +
 geom_point() +
 scale_x_continuous(name = "Proportion of genes DE") +
 scale_y_continuous(name = "-log10 P-value Pathways") +
 theme_bw()
dev.off()
cor(tab.path_prad$ DE_prop, -log10(tab.path_prad$P.Value), use = "complete")
# [1] 0.3691912

cor(tab.path_prad$DE_prop, abs(tab.path_prad$logFC), use = "complete")
# 0.3467721


cor(tab.path_prad$DE_prop[tab.path_prad$DE_prop > 0 ], abs(tab.path_prad$logFC)[tab.path_prad$DE_prop > 0 ], use = "complete")
# 0.345959
# comb.df <- left_join(tab.path, gos.mod, by = "category")
# comb.df$min_p <- pmin(comb.df$over_represented_pvalue, comb.df$under_represented_pvalue)
# plot(-log10(comb.df$min_p), -log10(comb.df$P.Value ))
# cor(-log10(comb.df$min_p), -log10(comb.df$P.Value ), use = "complete")

## GSVA
path_genes <- mclapply(paths.vec, function(x) subset(path.map, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec
prad_gsva <- gsva(prad, path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")


lm.gsva <- lmFit(assay(prad_gsva), mod) %>% eBayes()
tab.gsva_prad <- topTable(lm.gsva, coef = 2, n = Inf)
tab.gsva_prad$category <- rownames(tab.gsva_prad)

save(tab.gsva_prad, prad_gsva, file = "results/TCGA_PRAD/GSVA_results.Rdata")

prad_gsva.all <- gsva(prad.all, path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")
prad_gsva.all.filt <- prad_gsva.all[, colnames(prad)]

prad_gsva.mini <- gsva(prad[, prad$gleason == "High"], path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")

prad_gsva.control <- gsva(prad.all[, prad.all$sample_type == "Solid Tissue Normal"], path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")
cot <- cor(t(assay(prad_gsva.control)), t(assay(prad_gsva.all[, prad.all$sample_type == "Solid Tissue Normal"])))
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
hip.prad_vals <- get_paths_data(hip.res_prad )
hip.comp_prad <- do_wilcoxon(hip.prad_vals, hip.prad_vals$gleason, g1 = "High", g2 = "Low")

save(hip.comp_prad, hip.prad_vals, file = "results/TCGA_PRAD/hipathia.res.Rdata")


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
