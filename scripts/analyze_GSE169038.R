#'#################################################################################
#'#################################################################################
#' Analyze GSE169038 datasete
#'#################################################################################
#'#################################################################################


## Load libraries
library(fgsea)
library(limma)
library(tidyverse)
library(cowplot)
library(HDF5Array)
library(SummarizedExperiment)
library(e1071)
library(hipathia)
library(org.Hs.eg.db)
library(rjson)

load("results/GSE169038/allGenes.se.RData")
se.tcga_genes <- loadHDF5SummarizedExperiment("results/GSE169038/", prefix = "network_genes")


genes <- read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")
path.map <- read.table("results/preprocess/go_kegg_final_gene_map.tsv", header = TRUE)

prad.feat <- read.table("results/GSE169038/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
paths <- read.table("results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(prad.feat) <- paths.vec


pc.feat <- prcomp(prad.feat)

se$race <- ifelse(grepl( "White", se$characteristics_ch1.4), "EUR", "AFR")
se$decipher <- factor(gsub("decipher risk group: ", "", se$characteristics_ch1.3), levels = c("Lower", "Average", "Higher"))
se$primary <- gsub("primary gleason: ", "", se$characteristics_ch1.1)
se$secondary <- gsub("secondary gleason: ", "", se$characteristics_ch1.2)
se$primary <- as.numeric(ifelse(se$primary == "--", 1, se$primary))
se$secondary <- as.numeric(ifelse(se$secondary == "--", 1, se$secondary))
se$gleason_cat <- paste(se$primary, se$secondary, sep = "-")
se$gleason <- ifelse(se$primary == 5 |  se$secondary == 5 | se$gleason_cat == "4+4", "High", "Low")

## Subset samples with gleason < 3
se.filt <- se[, !(se$primary == 1 |  se$secondary == 1)]
prad.feat.filt <- prad.feat[!(se$primary == 1 |  se$secondary == 1), ]


## DE genes
mod <- model.matrix(~  gleason + race + decipher, colData(se.filt))
lm.genes <- lmFit(assay(se.filt), mod) %>% eBayes()
tab.genes_geo <- topTable(lm.genes, coef = 2, n = Inf)

tab.genes_geo$gene <- rowData(se)[as.character(rownames(tab.genes_geo )), "gene"]

save(tab.genes_geo, file = "results/GSE169038/de_genes_results.Rdata")


## fgsea
ranks <- tab.genes_geo$logFC
# names(ranks) <- rowData(se)[rownames(tab.genes), "gene"]
#
# names(paths.vec) <- paths.vec
# pathways <- lapply( paths.vec, function(x) subset(path.map, PathwayID  == x )$Symbol)
# fgseaRes <- fgsea(pathways = pathways, stats = ranks)
# fgseaRes$category <- fgseaRes$pathway

entrez_ids <- mapIds(org.Hs.eg.db, tab.genes_geo$gene, keytype="ENSEMBL", column="ENTREZID")

names(ranks) <- entrez_ids
pathways <- reactomePathways(as.character(entrez_ids))
# names(paths.vec) <- paths.vec
# pathways <- lapply( paths.vec, function(x) subset(path.map, PathwayID  == x )$Symbol)
fgseaRes_geo <- fgsea(pathways = pathways, stats = ranks)
save(fgseaRes_geo, file = "results/GSE169038/fgsea_results.Rdata")

## DE paths
lm.paths <- lmFit(t(prad.feat.filt), mod) %>% eBayes()
tab.paths_geo <- topTable(lm.paths, coef = 2, n = Inf)
tab.paths_geo$category <- rownames(tab.paths)
# tab.paths$pathway <- rownames(tab.paths)

# comb_fgsea <- left_join(fgseaRes, tab.paths, by = "pathway")

png("figures/GSE169038_pval_gseavsPaths.png")

ggplot(comb_fgsea, aes(x = -log10(P.Value ), y = -log10(pval ))) +
 geom_point() +
 scale_y_continuous(name = "-log10 P-value GSEA") +
 scale_x_continuous(name = "-log10 P-value Pathways") +
 theme_bw()
dev.off()

tab.paths_geo$DE_prop <- sapply( tab.paths_geo$category, function(cat) {
  genes <- subset(path.map, PathwayID  == cat )$Symbol
  mini_tab <- subset(tab.genes_geo, gene %in% genes)
  mean(mini_tab$adj.P.Val < 0.05)
})

png("figures/GSE169038_propDE_vs_pvalPaths.png")
ggplot(tab.paths_geo, aes(x = DE_prop, y = -log10(P.Value ))) +
 geom_point() +
 scale_x_continuous(name = "Proportion of genes DE") +
 scale_y_continuous(name = "-log10 P-value Pathways") +
 theme_bw()
dev.off()

png("figures/GSE169038_propDE_vs_logFCPaths.png")
ggplot(tab.paths_geo, aes(x = DE_prop , y = abs(logFC))) +
 geom_point() +
 scale_x_continuous(name = "Proportion of genes DE") +
 scale_y_continuous(name = "logFC Pathways (absolute value)") +
 theme_bw()
dev.off()


cor(tab.paths_geo$DE_prop, -log10(tab.paths_geo$P.Value))
# 0.2668079

cor(abs(tab.paths_geo$logFC), tab.paths_geo$DE_prop)
# [1] 0.2324701

cor(tab.paths_geo$DE_prop[tab.paths_geo$DE_prop > 0 ], abs(tab.paths_geo$logFC)[tab.paths_geo$DE_prop > 0 ], use = "complete")
# 0.229733



save(tab.paths_geo, file = "results/GSE169038/pathways_results.Rdata")



## hipathia
rownames(se) <- rowData(se)$gene
trans_data <- translate_data(se, "hsa")
exp_data <- normalize_data(trans_data)
hip_pathways <- load_pathways(species = "hsa")

hip.res_geo <- hipathia(exp_data, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.path_vals <- get_paths_data(hip.res_geo )
hip.comp_geo <- do_wilcoxon(hip.path_vals, hip.path_vals$gleason, g1 = "High", g2 = "Low")

save(hip.comp_geo, file = "results/GSE169038/hipathia.res.Rdata")



# Comparison between PRAD and GEO
## Pathways
load("results/TCGA_PRAD/pathways_results.Rdata")
comb_paths <- left_join(tab.path_prad, tab.paths_geo, by = "category", suffix = c(".TCGA", ".GEO")) %>%
  mutate(Signif = ifelse(adj.P.Val.TCGA < 0.05, ifelse(adj.P.Val.GEO < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val.GEO < 0.05, "GEO", "None")))
path.plot <- ggplot(comb_paths, aes(x = logFC.TCGA, y = logFC.GEO, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("GO + KEGG model") +
  xlab("logFC in TCGA") +
  ylab("logFC in GEO") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb_paths), cor(comb_paths$logFC.TCGA, comb_paths$logFC.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))



png("figures/TCGAvsGEO_logFC.png")
path.plot
dev.off()

png("figures/TCGAvsGEO_pval.png")
ggplot(comb_paths, aes(x = -log10(P.Value.TCGA), y =  -log10(P.Value.GEO), col = Signif)) +
  geom_point() +
  theme_bw()
dev.off()


## Add GO names
term <- read.table("data/GO_terms/term.txt", sep = "\t", quote = "", comment.char = "", as.is = TRUE)
term <- term[, c(2, 4)]
colnames(term) <- c("Name", "category")



kegg.annot <- fromJSON(file = "data/kegg_pathways.json")
kegg.df <- lapply(kegg.annot$children, function(x) {
  top_cat <- x$name
  paths.df <- lapply(x$children, function(y){
    cat2 <- y$name
    paths <- sapply(y$children, function(z) z$name)
    data.frame(top_cat = top_cat, category = cat2, path = paths)
  })
  Reduce(rbind, paths.df) %>%
    as_tibble()
}) %>%
 Reduce(rbind, .)
kegg.df <- kegg.df %>%
  mutate(category = paste0("path:hsa", substring(path, 0, 5)),
          Name = gsub("^[0-9]*  ", "", path))

path_map <- rbind(term, dplyr::select(kegg.df, Name, category))
comb_paths_annot <- left_join(comb_paths, path_map, by = "category") %>%
  dplyr::select(category, Name, ends_with("TCGA"), ends_with("GEO"), Signif) %>%
  dplyr::select(-starts_with("AveExpr"), -starts_with("t"), -starts_with("B"))
write.table(comb_paths_annot, file = "results/GSE169038/paths_results_comb.txt",
  sep = "\t", quote = FALSE, row.names = FALSE )


summary(lm(logFC.TCGA ~ logFC.GEO, comb_paths ))
summary(lm(logFC.TCGA ~ logFC.GEO, comb_paths, subset = Signif != "None" ))

cor(comb_paths$logFC.TCGA, comb_paths$logFC.GEO)
# [1] 0.4954506
cor(comb_paths[comb_paths$Signif != "None", ]$logFC.TCGA, comb_paths[comb_paths$Signif != "None", ]$logFC.GEO)
# [1] 0.6059908
table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif)
prop.table(table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif), margin = 2)
prop.table(table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif != "None"), margin = 2)

png("figures/GSE169038_TCGA_propDE_vs_logFCPaths.png", height = 300)
rbind(mutate(tab.path_prad, Dataset = "TCGA"),
    mutate(tab.paths_geo, Dataset = "GEO")) %>%
    mutate(Dataset = factor(Dataset, levels = c("TCGA", "GEO"))) %>%
    ggplot(aes(x = DE_prop , y = abs(logFC))) +
    geom_point() +
    scale_x_continuous(name = "Proportion of genes DE") +
    scale_y_continuous(name = "logFC Pathways (absolute value)") +
    facet_wrap(~ Dataset) +
    theme_bw()
dev.off()



## Compare Gene DE
load("results/TCGA_PRAD/genes_results.Rdata")
res_prad$gene <- rownames(res_prad)
comb.genes <- left_join(data.frame(res_prad), tab.genes_geo , by = "gene")   %>%
  mutate(Signif = ifelse(!is.na(padj) & padj  < 0.05, ifelse(adj.P.Val < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val < 0.05, "GEO", "None"))) %>%
  filter(!is.na(pvalue ) & !is.na(P.Value  ))

gene.plot <- ggplot(comb.genes, aes(x = log2FoldChange, y = logFC, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("Genes") +
  xlab("log2FC in TCGA") +
  ylab("logFC in GEO") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.genes), cor(comb.genes$log2FoldChange, comb.genes$logFC)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))

png("figures/TCGAvsGEO_genes_logFC.png")
gene.plot
dev.off()


cor(comb.genes$log2FoldChange, comb.genes$logFC)
# [1] 0.2488376
cor(comb.genes[comb.genes$Signif != "None", ]$log2FoldChange, comb.genes[comb.genes$Signif != "None", ]$logFC, use = "complete")
# [1] 0.3244926

table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif)
prop.table(table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif), margin = 2)
prop.table(table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif != "None"), margin = 2)


## GSEA
load("results/GSE169038/fgsea_results.Rdata")
comb.fgsea <- inner_join(fgseaRes_prad, fgseaRes_geo, by = "pathway", suffix = c(".TCGA", ".GEO"))   %>%
  mutate(Signif = ifelse(!is.na(padj.TCGA) & padj.TCGA  < 0.05, ifelse(padj.GEO < 0.05, "Both", "TCGA"),
                              ifelse(padj.GEO < 0.05, "GEO", "None")))


#
fgsea.plot <- ggplot(comb.fgsea, aes(x = NES.TCGA, y = NES.GEO, col = Signif)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("GSEA") +
  xlab("NES in TCGA") +
  ylab("NES in GEO") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.fgsea), cor(comb.fgsea$NES.TCGA, comb.fgsea$NES.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))

png("figures/TCGAvsGEO_GSEA_NES.png")
fgsea.plot
dev.off()

cor(comb.fgsea$NES.TCGA, comb.fgsea$NES.GEO)
# [1] 0.3988817
cor(comb.fgsea[comb.fgsea$Signif != "None", ]$NES.TCGA, comb.fgsea[comb.fgsea$Signif != "None", ]$NES.GEO, use = "complete")
# [1] 0.5000755

table(sign(comb.fgsea$NES.TCGA) == sign(comb.fgsea$NES.GEO), comb.fgsea$Signif)
prop.table(table(sign(comb.fgsea$NES.TCGA) == sign(comb.fgsea$NES.GEO), comb.fgsea$Signif), margin = 2)
prop.table(table(sign(comb.fgsea$NES.TCGA) == sign(comb.fgsea$NES.GEO), comb.fgsea$Signif != "None"), margin = 2)

## hipathia
load("results/TCGA_PRAD/hipathia.res.Rdata")
comb.hipathia <- inner_join(hip.comp_prad, hip.comp_geo, by = "name", suffix = c(".TCGA", ".GEO"))   %>%
  mutate(Signif = ifelse(!is.na(FDRp.value.TCGA) & FDRp.value.TCGA  < 0.05, ifelse(FDRp.value.GEO < 0.05, "Both", "TCGA"),
                              ifelse(FDRp.value.GEO < 0.05, "GEO", "None")))


#
hip.plot <- ggplot(comb.hipathia, aes(x = statistic.TCGA, y = statistic.GEO, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(name = "Significance", values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("Hipathia") +
  xlab("U statistic in TCGA") +
  ylab("U statistic in GEO") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.hipathia), cor(comb.hipathia$statistic.TCGA, comb.hipathia$statistic.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 20))

png("figures/TCGAvsGEO_hipathia_stat.png")
hip.plot
dev.off()

cor(comb.hipathia$statistic.TCGA, comb.hipathia$statistic.GEO)
# [1] 0.234694
cor(comb.hipathia[comb.hipathia$Signif != "None", ]$statistic.TCGA, comb.hipathia[comb.hipathia$Signif != "None", ]$statistic.GEO, use = "complete")
# [1] 0.4226619

table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif)
prop.table(table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif), margin = 2)
prop.table(table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif != "None"), margin = 2)

legend <- get_legend(
  # create some space to the left of the legend
  hip.plot + theme(legend.box.margin = margin(0, 0, 0, 12),
                    text = element_text(size = 25))+
                    guides(color = guide_legend(override.aes = list(size = 8)))
)


png("figures/TCGAvsGEO_panel.png", width = 950, height = 775)
plot_grid(
  plot_grid(gene.plot, fgsea.plot, hip.plot + theme(legend.position = "none"), path.plot, ncol = 2, labels = LETTERS[1:4], label_size = 20),
  legend, ncol = 2, rel_widths = c(5, 1))
dev.off()


## Compare values between TCGA and GEO
tcga.feat.all <- read.table("results/TCGA_gexp_coding_PRAD/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(tcga.feat.all) <- paths.vec

comb.feat.all <- rbind(tcga.feat.all, prad.feat)
pc.comb.all <- prcomp(comb.feat.all)

## Subset data
load("data/tcga_gexp_combat.Rdata")
tcga.prad <- gexp_tcga_combat[as.character(genes$V1), gexp_tcga_combat$project_id == "TCGA-PRAD"]


comb.all.df.pc <- data.frame(pc.comb.all$x, dataset = rep(c("TCGA", "GEO"), c(nrow(tcga.feat.all), nrow(prad.feat))),
  Type = c(tcga.prad$sample_type, rep("Primary Tumor", nrow(prad.feat) ))) %>%
  mutate(category = paste(dataset, Type))

png("figures/TCGA_GEO_all_samples_PCA.png")
ggplot(comb.all.df.pc, aes(x = PC1, y = PC2, color = category)) +
  geom_point() +
  theme_bw()
dev.off()

pc.tcga <- prcomp(tcga.feat)
tcga.df.pc <- data.frame(pc.tcga$x,  Type = tcga.prad$sample_type)
ggplot(tcga.df.pc, aes(x = PC1, y = PC2, color = Type)) +
  geom_point() +
  theme_bw()


tcga.feat <- read.table("results/TCGA_gexp_coding_PRAD_tumor/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(tcga.feat) <- paths.vec

tcga.prad <- tcga.prad[, !is.na(tcga.prad$paper_Reviewed_Gleason_category)]

## define TCGA gleason
tcga.prad$gleason <- ifelse(tcga.prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low")

comb.feat <- rbind(tcga.feat, prad.feat.filt)
comb.feat$dataset <- rep(c("TCGA", "GEO"), c(nrow(tcga.feat), nrow(prad.feat.filt)))
comb.feat$gleason <- c(tcga.prad$gleason, se.filt$gleason)


png("figures/TCGAvsGEO_GO0000212_boxplot.png")
ggplot(comb.feat, aes(x = gleason, y = `GO:0000212`, color = dataset)) +
  geom_boxplot() +
  theme_bw()
dev.off()




pc.comb <- prcomp(data.matrix(comb.feat[, 1:1337]))

comb.df.pc <- data.frame(pc.comb$x, dataset = rep(c("TCGA", "GEO"), c(nrow(tcga.feat), nrow(prad.feat.filt))),
      gleason = c(tcga.prad$gleason, se.filt$gleason))

png("figures/TCGAGEO_comb_PC_dataset.png")
ggplot(comb.df.pc, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point() +
  theme_bw()
dev.off()
#
png("figures/TCGAGEO_comb_PC_gleason.png")
ggplot(comb.df.pc, aes(x = PC1, y = PC2, color = gleason)) +
  geom_point() +
  theme_bw()
dev.off()



## Test SVM to classify gleason
load("results/TCGA_PRAD/svm_model.Rdata")
colnames(prad.feat.filt) <- gsub(":", "_", colnames( prad.feat.filt))
pred.geo <- predict(svm_gleason, prad.feat.filt)
table(prediction = pred.geo , real = se.filt$gleason )

pred.geo_filt <- predict(svm_gleason_filt, prad.feat.filt)
table(prediction = pred.geo_filt , real = se.filt$gleason )


## Train SVM to classify gleason
### All features
df_svm <-  data.frame(gleason = factor(se.filt$gleason), prad.feat.filt)
svm_gleason_geo <- svm(gleason ~ ., df_svm)
pred.geo2 <- predict(svm_gleason_geo, prad.feat.filt)
table(prediction = pred.geo2 , real = se.filt$gleason )

## Selected features
sel_paths <- subset(tab.paths, adj.P.Val < 0.05)
svm_gleason_geo_filt <- svm(gleason ~ ., df_svm[, c("gleason", gsub(":", "_", rownames( sel_paths)))])
pred.geo2_filt <- predict(svm_gleason_geo_filt, prad.feat.filt)
table(prediction = pred.geo2_filt , real = se.filt$gleason )

save(svm_gleason_geo, svm_gleason_geo_filt, file = "results/TCGA_PRAD/svm_model_geo.Rdata")
