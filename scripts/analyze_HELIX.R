#'#################################################################################
#'#################################################################################
#' Analyze HELIX datasete
#'#################################################################################
#'#################################################################################


## Load libraries
library(fgsea)
library(limma)
library(tidyverse)
library(cowplot)
library(HDF5Array)
library(SummarizedExperiment)
library(omicRexposome)
library(MultiDataSet)
library(rexposome)
library(mediation)
library(pheatmap)

load("results/HELIX/allGenes.se.RData")
se.tcga_genes <- loadHDF5SummarizedExperiment("results/HELIX/", prefix = "network_genes")


genes <- read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")
path.map <- read.table("results/preprocess/go_kegg_final_gene_map.tsv", header = TRUE)

helix.feat <- read.table("results/HELIX/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
paths <- read.table("results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(helix.feat) <- paths.vec

load("data/HELIX/exposome.RData")
description <- read_delim("data/HELIX/codebook_chf.csv", delim = ";" ) %>% data.frame()
rownames(description) <- description$variable_name

rownames(exposome) <- exposome$ID
exposome$ID <- NULL

## Merge phenotypes
rownames(phenotype) <- as.character(phenotype$ID)
colData(se.tcga_genes) <- cbind(colData(se.tcga_genes),  phenotype[colnames(se.tcga_genes), ])

## Test pathways for all the phenotypes
mod.asthma <- model.matrix(~ hs_asthma + e3_sex + age_sample_years + h_ethnicity_cauc, colData(se.tcga_genes))
lm.paths.asthma <- lmFit(t(helix.feat), mod.asthma ) %>% eBayes()
tab.paths.asthma <- topTable(lm.paths.asthma, coef = 2, n = Inf)
tab.paths.asthma$category <- rownames(tab.paths.asthma)

mod.asthma_cauc <- model.matrix(~ hs_asthma + e3_sex + age_sample_years, colData(se.tcga_genes[, se.tcga_genes$h_ethnicity_cauc == "yes"]))
lm.paths.asthma_cauc <- lmFit(t(helix.feat)[, se.tcga_genes$h_ethnicity_cauc == "yes"],mod.asthma_cauc  ) %>% eBayes()
tab.paths.asthma_cauc <- topTable(lm.paths.asthma_cauc, coef = 2, n = Inf)
tab.paths.asthma_cauc$category <- rownames(tab.paths.asthma_cauc)


mod.bw <- model.matrix(~ e3_bw + e3_sex + age_sample_years + h_ethnicity_cauc, colData(se.tcga_genes))
lm.paths.bw <- lmFit(t(helix.feat), mod.bw) %>% eBayes()
tab.paths.bw <- topTable(lm.paths.bw, coef = 2, n = Inf)
tab.paths.bw$category <- rownames(tab.paths.bw)

mod.bw_cauc <- model.matrix(~ e3_bw + e3_sex + age_sample_years, colData(se.tcga_genes[, se.tcga_genes$h_ethnicity_cauc == "yes"]))
lm.paths.bw_cauc <- lmFit(t(helix.feat)[, se.tcga_genes$h_ethnicity_cauc == "yes"],mod.bw_cauc  ) %>% eBayes()
tab.paths.bw_cauc <- topTable(lm.paths.bw_cauc, coef = 2, n = Inf)
tab.paths.bw_cauc$category <- rownames(tab.paths.bw_cauc)



mod.iq <- model.matrix(~ hs_correct_raven + e3_sex + age_sample_years + h_ethnicity_cauc, colData(se.tcga_genes))
lm.paths.iq <- lmFit(t(helix.feat), mod.iq ) %>% eBayes()
tab.paths.iq <- topTable(lm.paths.iq, coef = 2, n = Inf)
tab.paths.iq$category <- rownames(tab.paths.iq)

mod.neuro <- model.matrix(~ hs_Gen_Tot + e3_sex + age_sample_years + h_ethnicity_cauc, colData(se.tcga_genes))
lm.paths.neuro <- lmFit(t(helix.feat), mod.neuro ) %>% eBayes()
tab.paths.neuro <- topTable(lm.paths.neuro, coef = 2, n = Inf)
tab.paths.neuro$category <- rownames(tab.paths.neuro)

mod.bmi <- model.matrix(~ hs_zbmi_who + e3_sex + age_sample_years + h_ethnicity_cauc, colData(se.tcga_genes))
lm.paths.bmi <- lmFit(t(helix.feat), mod.bmi) %>% eBayes()
tab.paths.bmi <- topTable(lm.paths.bmi, coef = 2, n = Inf)
tab.paths.bmi$category <- rownames(tab.paths.bmi)
bmi_paths_sig <- subset(tab.paths.bmi, adj.P.Val < 0.05)$category

se.tcga_genes$obese <- ifelse(se.tcga_genes$hs_bmi_c_cat %in% c(1, 2), "normal", "obese")
mod.obese <- model.matrix(~ obese + e3_sex + age_sample_years + h_ethnicity_cauc, colData(se.tcga_genes))
lm.paths.obese <- lmFit(t(helix.feat), mod.obese) %>% eBayes()
tab.paths.obese <- topTable(lm.paths.obese, coef = 2, n = Inf)
tab.paths.obese$category <- rownames(tab.paths.obese)


## Test genes for all the phenotypes
lm.genes.asthma <- lmFit(assay(se), mod.asthma ) %>% eBayes()
tab.genes.asthma <- topTable(lm.genes.asthma, coef = 2, n = Inf)
tab.genes.asthma$category <- rownames(tab.genes.asthma)

lm.genes.bw <- lmFit(assay(se), mod.bw ) %>% eBayes()
tab.genes.bw <- topTable(lm.genes.bw, coef = 2, n = Inf)
tab.genes.bw$category <- rownames(tab.genes.bw)

lm.genes.bw_cauc <- lmFit(assay(se[, se$h_ethnicity_cauc == "yes"]), mod.bw_cauc ) %>% eBayes()
tab.genes.bw_cauc <- topTable(lm.genes.bw_cauc, coef = 2, n = Inf)
tab.genes.bw_cauc$gene <- rownames(tab.genes.bw)

lm.genes.iq <- lmFit(assay(se), mod.iq ) %>% eBayes()
tab.genes.iq <- topTable(lm.genes.iq, coef = 2, n = Inf)
tab.genes.iq$category <- rownames(tab.genes.iq)

lm.genes.neuro <- lmFit(assay(se), mod.neuro ) %>% eBayes()
tab.genes.neuro <- topTable(lm.genes.neuro, coef = 2, n = Inf)
tab.genes.neuro$category <- rownames(tab.genes.neuro)

lm.genes.bmi <- lmFit(assay(se), mod.bmi ) %>% eBayes()
tab.genes.bmi <- topTable(lm.genes.bmi, coef = 2, n = Inf)
tab.genes.bmi$category <- rownames(tab.genes.bmi)
bmi_genes_sig <- subset(tab.genes.bmi, adj.P.Val < 0.05)$category

lm.genes.obese <- lmFit(assay(se), mod.obese ) %>% eBayes()
tab.genes.obese <- topTable(lm.genes.obese, coef = 2, n = Inf)
tab.genes.obese$category <- rownames(tab.genes.obese)


## Compute PCs
pc.feat <- prcomp(helix.feat)
pc.genes <- prcomp(t(assay(se)))


### Overlap between bmi and bw
png("figures/HELIX.bmi_bw.png", width = 1000)
par(mfrow = c(1, 2))
plot(tab.paths.bmi$logFC, tab.paths.bw[rownames(tab.paths.bmi), "logFC"], main ="Pathways")
plot(tab.genes.bmi$logFC, tab.genes.bw[rownames(tab.genes.bmi), "logFC"], main ="Genes")
dev.off()

# Association with exposures

## ExWAS
pheno_comb <- left_join(phenotype %>% mutate(ID = as.character(ID)), data.frame(colData(se)) , by = "ID")
expo_set <- loadExposome(exposome, description[colnames(exposome), ], pheno_comb)
expo_set_post <- expo_set[fData(expo_set)$period == "Postnatal",]

bmi_expos <- exwas(expo_set_post, formula = hs_zbmi_who ~ e3_sex + age_sample_years + h_ethnicity_cauc, family = "gaussian")
bmi_expos_df <- extract(bmi_expos)
bmi_expos_sig <- rownames(subset(bmi_expos_df, pvalue < tef(bmi_expos) ))


# bw_expos <- exwas(expo_set_post, formula = e3_bw ~ e3_sex + age_sample_years + h_ethnicity_cauc, family = "gaussian")
# bw_expos_df <- extract(bw_expos)
# bw_expos_sig <- rownames(subset(bw_expos_df, pvalue < tef(bw_expos) ))

## Create MDS
eset_genes <- ExpressionSet(assay(se), AnnotatedDataFrame(data.frame(colData(se))),
             AnnotatedDataFrame(data.frame(rowData(se))))
colnames(fData(eset_genes))[colnames(fData(eset_genes)) == "seqname"] <- "chromosome"
colnames(fData(eset_genes))[colnames(fData(eset_genes)) == "stop"] <- "end"

mds_genes <- createMultiDataSet()
mds_genes <- add_genexp(mds_genes, eset_genes)
mds_genes <- add_exp(mds_genes, expo_set_post)

rownames(helix.feat) <- colnames(se)

fdata <- data.frame(chromosome = rep("1", 1337), start = rep(1, 1337), end = rep(1, 1337))
rownames(fdata) <- colnames(helix.feat)

eset_paths <- ExpressionSet(t(helix.feat), AnnotatedDataFrame(data.frame(colData(se))),
  AnnotatedDataFrame(fdata))


mds_paths <- createMultiDataSet()
mds_paths <- add_genexp(mds_paths, eset_paths)
mds_paths <- add_exp(mds_paths, expo_set_post)

## Association with genes
assoc_gene <- association(mds_genes, formula = ~ e3_sex + age_sample_years + h_ethnicity_cauc,
    expset = "exposures", omicset = "expression")

save(assoc_gene, file = "results/HELIX/genes_expos_assoc.Rdata")

#
hit_gene <- tableHits(assoc_gene, th = 0.05/nrow(se))
lambda_gene <- tableLambda(assoc_gene)
merge(hit_gene, lambda_gene, by="exposure")


assoc_gene_sig <- Reduce(rbind, lapply(names(assoc_gene), function(expo){
  mini <- getAssociation(assoc_gene, rid = expo)
  mini$exposure <- expo
  mini$gene <- rownames(mini)

  filter(mini, adj.P.Val < 0.05)
  })
)

## Association with pathways
assoc_path <- association(mds_paths, formula = ~ e3_sex + age_sample_years + h_ethnicity_cauc,
    expset = "exposures", omicset = "expression")
save(assoc_path, file = "results/HELIX/path_expos_assoc.Rdata")
#
hit_path <- tableHits(assoc_path, th = 0.05/nrow(eset_paths))
lambda_path <- tableLambda(assoc_path)
merge(hit_path, lambda_path, by="exposure")


hit_comb <- left_join(hit_path, hit_gene, by = "exposure", suffix = c(".path", ".gene")) %>%
    left_join(dplyr::select(description, variable_name, family, domain), by = c("exposure" = "variable_name"))

png("figures/HELIX_exwas_genes_paths.png", width = 800)
ggplot(hit_comb, aes(x = hits.path, y = hits.gene, color = family)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = nrow(se)/nrow(eset_paths)) +
  facet_wrap(~ domain) +
  xlab("DE pathways (Bonferroni)") +
  ylab("DE genes (Bonferroni)")
dev.off()

cor(hit_comb$hits.path, hit_comb$hits.gene)
# [1] 0.8232468

assoc_path_sig <- Reduce(rbind, lapply(names(assoc_path), function(expo){
  mini <- getAssociation(assoc_path, rid = expo)
  mini$exposure <- expo
  mini$path <- rownames(mini)

  filter(mini, adj.P.Val < 0.05)
  })
)

## Search association in bmi
subset(assoc_path_sig, exposure %in% bmi_expos_sig & path %in% bmi_paths_sig)
subset(assoc_gene_sig, exposure %in% bmi_expos_sig & gene %in% bmi_genes_sig)

df.med_bmi <- data.frame(expo = unlist(expos(expo_set_post["hs_cu_c_Log2", colnames(eset_paths)])),
                 path = as.numeric(exprs(eset_paths["GO:0032688", ] )))
df.med_bmi <- cbind(df.med_bmi, pData(eset_paths)) %>%
  left_join(phenotype %>% mutate(ID = as.character(ID)), by = "ID")


mod.med <- lm(path ~  expo + e3_sex + age_sample_years + h_ethnicity_cauc, df.med_bmi )
mod.out <- lm(hs_zbmi_who  ~ path + expo + e3_sex + age_sample_years + h_ethnicity_cauc, df.med_bmi)
med <- mediate(mod.med, mod.out, treat = "expo", mediator = "path", sims = 1000)


expo_path_plot <- ggplot(df.med_bmi, aes(x = expo, y = path)) +
  geom_point() +
  theme_bw() +
  xlab("Cu") +
  ylab("GO:0032688") +
  geom_smooth(method = "lm")



expo_bmi_plot <- ggplot(df.med_bmi, aes(x = expo, y = hs_zbmi_who)) +
  geom_point() +
  theme_bw() +
  xlab("Cu") +
  ylab("zBMI") +
  geom_smooth(method = "lm")

path_bmi_plot <- ggplot(df.med_bmi, aes(x = path, y = hs_zbmi_who)) +
  geom_point() +
  theme_bw() +
  xlab("GO:0032688") +
  ylab("zBMI") +
  geom_smooth(method = "lm")


png("figures/HELIX_cu_bmi_go.png", width = 1200)
plot_grid(expo_path_plot, path_bmi_plot, expo_bmi_plot, labels = LETTERS[1:3], nrow = 1)
dev.off()

## Association with genes from pathway
path_genes <- subset(path.map, PathwayID == "GO:0032688")$Symbol

path_tcs <- rowData(se.tcga_genes[path_genes[path_genes %in% rownames(se.tcga_genes)], ])$transcript_id

relevance <- h5read("results/HELIX/comb_paths3_v3.6/relevance/relevance_genes_GO0032688.h5","GO0032688")
genes <- as.character(read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")$V1)
rownames(relevance) <- genes
path_genesf <- path_genes[path_genes %in% genes]

weights <- h5read("results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/model_weights.h5","weights_paths")
rownames(weights) <- paths.vec
colnames(weights) <- genes


relevance_path <- relevance[path_genesf, ]
pc_rel <- prcomp(relevance_path)

runMediations <- function(gene){
  df.med_bmi$tc <- as.numeric(assay(se[as.character(gene), ]))
  mod.med <- lm(tc ~  expo + e3_sex + age_sample_years + h_ethnicity_cauc, df.med_bmi )
  mod.out <- lm(hs_zbmi_who  ~ tc + expo + e3_sex + age_sample_years + h_ethnicity_cauc, df.med_bmi)
  med <- mediate(mod.med, mod.out, treat = "expo", mediator = "tc", sims = 1000)
  med
}


med_genes_path <- lapply(path_tcs, runMediations)


cor(t(assay(se[cu_genes,])), df.med_bmi$path)
cor(t(assay(se[path_tcs,])), df.med_bmi$path)

cu_assoc_genes <- getAssociation(assoc_gene, rid = "hs_cu_c_Log2")
path_genes_df <- data.frame(ensembl = rowData(se[path_tcs, ])$gene,
                            r = cor(t(assay(se[path_tcs,])), df.med_bmi$path),
                            weight = weights["GO:0032688", path_genesf],
                            logFC_expo = cu_assoc_genes[path_tcs, "logFC"],
                            pval_expo = cu_assoc_genes[path_tcs, "P.Value"],
                            adjpval_expo = cu_assoc_genes[path_tcs, "adj.P.Val"],
                            logFC_bmi = tab.genes.bmi[path_tcs, "logFC"],
                            pval_bmi = tab.genes.bmi[path_tcs, "P.Value"],
                            adjpval_bmi = tab.genes.bmi[path_tcs, "adj.P.Val"],
                            prop_med = sapply(med_genes_path, function(x) x$n.avg),
                            pval_med = sapply(med_genes_path, function(x) x$n.avg.p)) %>%
                mutate(Significance = ifelse(adjpval_expo < 0.05, ifelse(adjpval_bmi < 0.05, "Both", "Cu"),
                                                                  ifelse(adjpval_bmi < 0.05, "BMI", "None")),
                      Importance = ifelse(abs(r) > 0.2, "High importance", "Low importance")) %>%
                rbind(data.frame(ensembl = "GO:0032688", r = 1, weight = NA,
                                logFC_expo = assoc_path_sig[ "GO:0032688", "logFC"],
                                pval_expo = assoc_path_sig[ "GO:0032688", "P.Value"],
                                adjpval_expo = assoc_path_sig[ "GO:0032688", "adj.P.Val"],
                                logFC_bmi = tab.paths.bmi["GO:0032688", "logFC"],
                                pval_bmi = tab.paths.bmi["GO:0032688", "P.Value"],
                                adjpval_bmi = tab.paths.bmi["GO:0032688", "adj.P.Val"],
                                prop_med = med$n.avg,
                                pval_med = med$n.avg.p,
                              Significance = "Both", Importance = "GO:0032688"))


png("figures/HELIX_GO0032688_genes.png", width = 480, height = 400)
ggplot(path_genes_df, aes(x = logFC_bmi, y = logFC_expo, shape = Significance, color = Importance)) +
  geom_point() +
  theme_bw() +
  xlab("logFC for BMI") +
  ylab("logFC for Cu")
dev.off()

col_colors <- list(
  hs_bmi_c_cat = c("1" = "gray", "2" = "black", "3" = "pink", "4" = "red"),
  e3_sex = c("female" = "purple", "male" = "lightblue")
)


png("figures/HELIX_GO0032688_heatmap.png", width = 1000)
pheatmap(assay(se[rownames(subset(path_genes_df, Importance == "High importance")),]),  scale = "none",
         annotation_col  = pData(expo_set_post[, colnames(se)])[, c("hs_bmi_c_cat", "e3_sex")],
         annotation_colors =  col_colors,
         show_colnames = FALSE)
dev.off()

clus <- hclust(dist(t(assay(se[rownames(subset(path_genes_df, Importance == "High importance")),]))))
table(cutree(clus, 3), expo_set_post[, colnames(se)]$hs_bmi_c_cat)
#     1   2   3   4
 # 1   2 350 117  67
 # 2   7 354  75  34
 # 3   0   1   0   0

library(e1071)
cot <- data.frame(helix.feat)
cot$bmi <- expo_set_post[, colnames(se)]$hs_bmi_c_cat
model_svm <- svm(bmi ~ ., cot[1:500, ])

pred <- predict(model_svm, cot[, -1338])
table(prediction = pred, real = expo_set_post[, colnames(se)]$hs_bmi_c_cat, seq_len(nrow(cot)) %in% seq_len(500))


## Summary of genes associated with Cu and BMI
cu_genes <- subset(assoc_gene_sig, exposure == "hs_cu_c_Log2" & gene %in% bmi_genes_sig)$gene
med_genes <- lapply(cu_genes, runMediations)

cu_genes_df <- data.frame(gene = rowData(se[cu_genes, ])$GeneSymbolDB2,
                            logFC_expo = cu_assoc_genes[cu_genes, "logFC"],
                            pval_expo = cu_assoc_genes[cu_genes, "P.Value"],
                            adjpval_expo = cu_assoc_genes[cu_genes, "adj.P.Val"],
                            logFC_bmi = tab.genes.bmi[cu_genes, "logFC"],
                            pval_bmi = tab.genes.bmi[cu_genes, "P.Value"],
                            adjpval_bmi = tab.genes.bmi[cu_genes, "adj.P.Val"],
                            prop_med = sapply(med_genes, function(x) x$n.avg),
                            pval_med = sapply(med_genes, function(x) x$n.avg.p)) %>%
                mutate(Significance = ifelse(adjpval_expo < 0.05, ifelse(adjpval_bmi < 0.05, "Both", "Cu"),
                                                                  ifelse(adjpval_bmi < 0.05, "BMI", "None")))
write.table(cu_genes_df, file = "results/HELIX/cu_genes_assoc.txt", sep = "\t", quote = FALSE)

## Multidimensional
path_mcca <- crossomics(mds_paths, method = "mcca")
path_mcia <- crossomics(mds_paths, method = "mcia")
