#'#################################################################################
#'#################################################################################
#' Explore whether NetActivity improves the separation between control and tumor
#'#################################################################################
#'#################################################################################


## Load libraries
library(SummarizedExperiment)
library(HDF5Array)
library(DESeq2)
library(NetActivity)
library(NetActivityData)
library(tidyverse)
library(cowplot)
library(parallel)
library(vegan)
library(ggcorrplot)

data(gtex_gokegg)
input_genes <- read.table("results/GTEx/input_genes.txt")

gtex.vst <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_")
assay(gtex.vst) <- data.matrix(assay(gtex.vst))
preproc <- prepareSummarizedExperiment(gtex.vst, "gtex_gokegg")
gtex_scores <- computeGeneSetScores(preproc, "gtex_gokegg")
rda_gtex <- rda(t(assay(gtex.vst)) ~ t(assay(gtex_scores)))
save(rda_gtex, file = "results/manuscript/mainModel_GTEx_explainability.Rdata")


# TCGA
## Load  data
load("data/tcga_gexp_combat.Rdata")
## Restrict analysis to input genes
tcga_input <- gexp_tcga_combat[as.character(input_genes$V1), ]

## Overall variability
counts <- assay(tcga_input)
mode(counts) <- "integer"
counts[is.na(counts)] <- max(counts, na.rm = TRUE)
deseq <- DESeqDataSetFromMatrix(countData = counts,
              colData = colData(gexp_tcga_combat),
              design = ~ 1)
vst <-  vst(deseq, blind=FALSE)

preproc <- prepareSummarizedExperiment(vst, "gtex_gokegg")
all_scores <- computeGeneSetScores(preproc, "gtex_gokegg")
rda_tcga <- rda(t(assay(vst)) ~ t(assay(all_scores)))
save(rda_tcga, file = "results/manuscript/mainModel_TCGA_explainability.Rdata")

## Select tumors
tumor_tab <- table(gexp_tcga_combat$project_id)
sel_tumors <- names(tumor_tab)
names(sel_tumors) <- sel_tumors

## Compute PCs
tumor_pcs <- mclapply(sel_tumors, function(tum){
  message(tum)
  set <- tcga_input[, tcga_input$project_id == tum]

  counts <- assay(set)
  mode(counts) <- "integer"
  counts[is.na(counts)] <- max(counts, na.rm = TRUE)
  deseq <- DESeqDataSetFromMatrix(countData = counts,
                      colData = colData(set),
                      design = ~ 1)
  vst <-  vst(deseq, blind=FALSE)
  pc_vst <- prcomp(t(assay(vst)), rank. = 10)

  preproc <- prepareSummarizedExperiment(vst, "gtex_gokegg")
  scores <- computeGeneSetScores(preproc, "gtex_gokegg")
  pc_scores <- prcomp(t(assay(scores)), rank. = 10)

  list(
    raw = vst,
    vst = pc_vst,
    scores = pc_scores)
}, mc.cores = 5)
rda_pc10 <- mclapply(tumor_pcs, function(y) rda(t(assay(y$raw)) ~ y$scores$x), mc.cores = 10)
var_prop10 <- sapply(rda_pc10, function(i) i$CCA$tot.chi/i$tot.chi)
pc_prop10 <- sapply(tumor_pcs, function(x) sum(x$vst$sdev[1:10]**2)/sum(x$vst$sdev**2))

# #
# tumor_cors <- lapply(tumor_pcs, function(l){
#   cor(l$vst$x, l$scores$x)
# })
# sapply(tumor_cors, function(x) x$vst[1,1])
#
# bad_cor <- names(which(abs(sapply(tumor_cors, function(x) x$vst[1,1])) < 0.7) )
# sapply(tumor_cors[bad_cor], function(x) x$vst_filt[1,1])

# rda_pc2 <- lapply(tumor_pcs, function(y) rda(y$vst$x[, 1:2] ~ y$scores$x[, 1:2]))


# df_pcs_prad <- data.frame(PC1 = c(tumor_pcs[["TCGA-PRAD"]]$vst$x[, 1], tumor_pcs[["TCGA-PRAD"]]$scores$x[, 1]),
#                      PC2 = c(tumor_pcs[["TCGA-PRAD"]]$vst$x[, 2], -tumor_pcs[["TCGA-PRAD"]]$scores$x[, 2]),
#                       dataset = rep(c("Original gene expression", "Gene set activities"), c(nrow(tumor_pcs[["TCGA-PRAD"]]$vst$x), nrow(tumor_pcs[["TCGA-PRAD"]]$scores$x))),
#                     Sample = rep(factor(tumor_pcs[["TCGA-PRAD"]]$raw$sample_type), 2)) %>%
#                     mutate(dataset = factor(dataset, levels = c("Original gene expression", "Gene set activities")),
#                           Sample = ifelse(Sample == "Solid Tissue Normal", "Normal", "Cancer"),
#                           Sample = factor(Sample, levels = c("Normal", "Cancer")))

# pcs_model <- ggplot(df_pcs_prad, aes(x = PC1, y = PC2, color = Sample)) +
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~ dataset, scales = "free") +
#   scale_color_manual(values = c("green", "red")) +
# ggtitle("PRAD dataset") +
#   theme(plot.title = element_text(hjust = 0.5))

#
pcs_model_prad_ori <- data.frame(tumor_pcs[["TCGA-PRAD"]]$vst$x,
                       Sample = factor(tumor_pcs[["TCGA-PRAD"]]$raw$sample_type)) %>%
                    mutate( Sample = ifelse(Sample == "Solid Tissue Normal", "Normal", "Cancer"),
                            Sample = factor(Sample, levels = c("Normal", "Cancer"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Sample)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = c("green", "red")) +
    ggtitle("PRAD - original gene expression") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    xlab(sprintf("PC1 (%.1f%%)", tumor_pcs[["TCGA-PRAD"]]$vst$sdev[1]**2/sum(tumor_pcs[["TCGA-PRAD"]]$vst$sdev**2)*100)) +
    ylab(sprintf("PC2 (%.1f%%)", tumor_pcs[["TCGA-PRAD"]]$vst$sdev[2]**2/sum(tumor_pcs[["TCGA-PRAD"]]$vst$sdev**2)*100))

#
pcs_model_prad_scores <- data.frame(tumor_pcs[["TCGA-PRAD"]]$scores$x,
                       Sample = factor(tumor_pcs[["TCGA-PRAD"]]$raw$sample_type)) %>%
                    mutate( Sample = ifelse(Sample == "Solid Tissue Normal", "Normal", "Cancer"),
                            Sample = factor(Sample, levels = c("Normal", "Cancer"))) %>%
  ggplot(aes(x = PC1, y = -PC2, color = Sample)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = c("green", "red")) +
    ggtitle("PRAD - NetActivity GSAS") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(sprintf("PC1 (%.1f%%)", tumor_pcs[["TCGA-PRAD"]]$scores$sdev[1]**2/sum(tumor_pcs[["TCGA-PRAD"]]$scores$sdev**2)*100)) +
    ylab(sprintf("PC2 (%.1f%%)", tumor_pcs[["TCGA-PRAD"]]$scores$sdev[2]**2/sum(tumor_pcs[["TCGA-PRAD"]]$scores$sdev**2)*100))

cor_prad <- cor(tumor_pcs[["TCGA-PRAD"]]$vst$x, tumor_pcs[["TCGA-PRAD"]]$scores$x)

corplot_prad <- ggcorrplot(abs(cor_prad), method = "circle", hc.order = FALSE,
      title = "PRAD") +
  scale_x_discrete("Original gene expression") +
  scale_y_discrete(name = "NetActivity GSAS") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
   axis.title.x = element_text(angle = 0, color = 'grey20', size = 14),
   axis.title.y = element_text(angle = 90, color = 'grey20', size = 14))



## IBD cohort
vst_ibd <- loadHDF5SummarizedExperiment("results/SRP042228/", prefix = "vsd_norm")
vst_ibd_input <- vst_ibd[as.character(input_genes$V1), ]
pc_vst_ibd <- prcomp(t(assay(vst_ibd_input)), rank. = 10)

assay(vst_ibd_input) <- data.matrix(assay(vst_ibd_input))
preproc_ibd <- prepareSummarizedExperiment(vst_ibd_input, "gtex_gokegg")
scores_ibd <- computeGeneSetScores(preproc_ibd, "gtex_gokegg")
pc_scores_ibd <- prcomp(t(assay(scores_ibd)), rank. = 10)
rda_ibd <- rda(t(assay(vst_ibd_input)) ~ pc_scores_ibd$x)

pc_summary <- data.frame(prop10 = c(var_prop10, rda_ibd$CCA$tot.chi/rda_ibd$tot.chi),
          dataset = c(names(var_prop10), "IBD"),
          pc10 = c(pc_prop10, sum(pc_vst_ibd$sdev[1:10]**2)/sum(pc_vst_ibd$sdev**2)))

# pc_plot1 <- ggplot(pc_summary, aes(x = abs(cor1)**2*100, y = pc1*100)) +
#     geom_point() +
#     theme_bw() +
#     ggrepel::geom_label_repel(aes(label = tumor)) +
#     xlab("Variance explained by scores (%)") +
#     ylab("Variance (%)") +
#     ggtitle("PC1 correlation") +
#     theme(plot.title = element_text(hjust = 0.5))
#
pc_plot10 <- ggplot(pc_summary, aes(x = prop10*100, y = pc10*100)) +
    geom_point() +
    theme_bw() +
    ggrepel::geom_label_repel(data = subset(pc_summary, dataset %in% c("TCGA-PRAD", "IBD")),
                              aes( label = dataset))  +
    xlab("Variance explained by 10 top GSAS (%)") +
    ylab("Variance explained by 10 top PCs (%)") +
      theme(plot.title = element_text(hjust = 0.5))

# png("figures/tumors_score_variance.png", width = 1300)
# plot_grid(pc_plot1, pc_plot10, ncol = 2)
# dev.off()

# png("figures/tumors_score_variance.png", height = 1200, width = 1200, res = 300)
# pc_plot10
# dev.off()


png("figures/tumors_score_variance_abstract.png", height = 300,  width = 600)
ggplot(pc_summary, aes(x = prop10*100, y = pc10*100)) +
    geom_point(size = 2.5) +
    theme_bw() +
    xlab("Variance explained by the PCs 1-10 of scores (%)") +
    ylab("Variance of PCs 1-10 (%)") +
    ggtitle("Scores representativity") +
      theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))
dev.off()

summary(lm(prop10 ~ pc10, pc_summary))


df_pcs_ibd <- data.frame(PC1 = c(pc_vst_ibd$x[, 1], pc_scores_ibd$x[, 1]),
                     PC2 = c(pc_vst_ibd$x[, 2], pc_scores_ibd$x[, 2]),
                      dataset = rep(c("Original gene expression", "GSAS"), c(nrow(pc_vst_ibd$x), nrow(pc_scores_ibd$x))),
                    Sample = rep(factor(vst_ibd$diagnosis), 2)) %>%
                    mutate(dataset = factor(dataset, levels = c("Original gene expression", "GSAS")),
                            Sample = recode(Sample, Control = "Normal"),
                            Sample = factor(Sample, levels = c("Normal", "CD", "UC")))



#
pcs_model_ibd <- ggplot(df_pcs_ibd, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_manual(values = c("green", "blue", "orange")) +
  ggtitle("IBD dataset") +
  theme(plot.title = element_text(hjust = 0.5))


#
pcs_model_ibd_ori <- data.frame(pc_vst_ibd$x,
                       Sample = factor(vst_ibd$diagnosis)) %>%
                    mutate( Sample = recode(Sample, Control = "Normal"),
                            Sample = factor(Sample, levels = c("Normal", "CD", "UC"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Sample)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = c("green", "blue", "orange")) +
    ggtitle("IBD - original gene expression") +
    theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none") +
    xlab(sprintf("PC1 (%.1f%%)", pc_vst_ibd$sdev[1]**2/sum(pc_vst_ibd$sdev**2)*100)) +
    ylab(sprintf("PC2 (%.1f%%)", pc_vst_ibd$sdev[2]**2/sum(pc_vst_ibd$sdev**2)*100))


#
pcs_model_ibd_scores <- data.frame(pc_scores_ibd$x,
                       Sample = factor(vst_ibd$diagnosis)) %>%
                    mutate( Sample = recode(Sample, Control = "Normal"),
                            Sample = factor(Sample, levels = c("Normal", "CD", "UC"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Sample)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = c("green", "blue", "orange")) +
    ggtitle("IBD - NetActivity GSAS") +
    theme(plot.title = element_text(hjust = 0.5))  +
    xlab(sprintf("PC1 (%.1f%%)", pc_scores_ibd$sdev[1]**2/sum(pc_scores_ibd$sdev**2)*100)) +
    ylab(sprintf("PC2 (%.1f%%)", pc_scores_ibd$sdev[2]**2/sum(pc_scores_ibd$sdev**2)*100))



corplot_ibd <- ggcorrplot(abs(cor(pc_vst_ibd$x, pc_scores_ibd$x)), method = "circle", hc.order = FALSE,
      title = "IBD") +
  scale_x_discrete("Original gene expression") +
  scale_y_discrete(name = "NetActivity GSAS") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
   axis.title.x = element_text(angle = 0, color = 'grey20', size = 14),
   axis.title.y = element_text(angle = 90, color = 'grey20', size = 14))

load("results/manuscript/mainModel_GTEx_explainability.Rdata")
load("results/manuscript/mainModel_TCGA_explainability.Rdata")

plot_represent <- data.frame(Perc = c(rda_gtex$CCA$tot.chi/rda_gtex$tot.chi, rda_tcga$CCA$tot.chi/rda_tcga$tot.chi), 
                             Dataset = c("GTEx", "TCGA")) %>%
  ggplot(aes(x = Dataset, y = Perc*100)) +
  geom_bar(stat = "identity") +
  xlab("Dataset") +
  ylab("Variance expl. by GSAS (%)") +
  theme_bw() +
  ylim(c(0, 100)) +
  theme(text = element_text(size = 14))


#
load("results/manuscript/df_gtex_training_replicability.Rdata")


plot_rep <- ggplot(df.path_sel, aes(x = group, y = minCor, color = training)) +
 geom_boxplot() +
 scale_x_discrete(name = "") +
 scale_y_continuous(name = "Replicability") +
 theme_bw() +
 scale_color_manual(name = "Training", values = c("grey", "#777777", "black")) +
 theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 14))

png("figures/evaluationTCGAPRAD_panel2.png", width = 4500, height = 3300, res = 300)
plot_grid(
  plot_grid(plot_rep, plot_represent, pc_plot10, rel_widths = c(7, 4, 4), nrow = 1, labels = c("A", "B", "C")),
  plot_grid(pcs_model_ibd_ori, pcs_model_ibd_scores, corplot_ibd, nrow = 1, rel_widths = c(5, 7, 5), labels = c("D", "E", "H")),
  plot_grid(pcs_model_prad_ori, pcs_model_prad_scores, corplot_prad, nrow = 1, rel_widths = c(5, 7, 5), labels = c("F", "G", "I")), 
  nrow = 3
)
dev.off()
