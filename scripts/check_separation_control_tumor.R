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

## Load TCGA data
load("data/tcga_gexp_combat.Rdata")
data(gtex_gokegg)

## Select tumors
tumor_tab <- table(gexp_tcga_combat$project_id)
sel_tumors <- names(tumor_tab)
names(sel_tumors) <- sel_tumors


## Compute PCs
tumor_pcs <- mclapply(sel_tumors, function(tum){
  message(tum)
  set <- gexp_tcga_combat[, gexp_tcga_combat$project_id == tum]

  counts <- assay(set)
  mode(counts) <- "integer"
  counts[is.na(counts)] <- max(counts, na.rm = TRUE)
  deseq <- DESeqDataSetFromMatrix(countData = counts,
                      colData = colData(set),
                      design = ~ 1)
  vst <-  vst(deseq, blind=FALSE)

  pc_vst <- prcomp(t(assay(vst)), rank. = 10)

  # pc_vst_filt <- prcomp(t(assay(vst[colnames(gtex_gokegg), ])), rank. = 10)


  preproc <- prepareSummarizedExperiment(vst, "gtex_gokegg")
  # pc_preproc <- prcomp(t(assay(preproc[colnames(gtex_gokegg), ])), rank. = 10)

  scores <- computeGeneSetScores(preproc, "gtex_gokegg")
  pc_scores <- prcomp(t(assay(scores)), rank. = 10)

  list(vst = pc_vst,
     # vst_filt = pc_vst_filt, preproc = pc_preproc,
     scores = pc_scores, status = set$sample_type)
}, mc.cores = 10)
# tumor_lms <- lapply(tumor_pcs, function(l){
#
#   stat <- ifelse( l$status == "Solid Tissue Normal", "Control", "Cancer")
#
#   if(all(stat == "Cancer")){
#     return(NULL)
#   }
#
#   list(vst = summary(lm(l$vst$x[, 1] ~ stat)),
#   vst_filt = summary(lm(l$vst_filt$x[, 1] ~ stat)),
#   preproc = summary(lm(l$preproc$x[, 1] ~ stat)),
#   scores = summary(lm(l$scores$x[, 1] ~ stat)))
#
# })
# cot <- sapply(tumor_lms, function(x) {
#   if (is.null(x)){
#     return(c(vst = NA, vst_filt = NA, preproc = NA, scores = NA))
#   }
#   c(vst = x$vst$r.squared, vst_filt = x$vst_filt$r.squared, preproc = x$preproc$r.squared, scores = x$scores$r.squared)
# }) %>%
#   as_tibble( )%>%
#   mutate(method = c("VST_all", "VST_model", "Standardized", "Scores")) %>%
#   gather(Tumor, R2, 1:33)
#
# ggplot(cot, aes(x = method, y = R2)) +
#   geom_boxplot() +
#   theme_bw()
#
# cot %>%
#     filter(!is.na(R2)) %>%
#     spread(method, R2) %>%
#     ggplot(aes(x = VST_model, y = Scores)) +
#     geom_point() +
#     theme_bw() +
#     geom_abline(slope = 1) +
#     xlim(c(0, 1)) +
#     ylim(c(0, 1))
#
#

#
tumor_cors <- lapply(tumor_pcs, function(l){
  cor(l$vst$x, l$scores$x)
})
# sapply(tumor_cors, function(x) x$vst[1,1])
#
# bad_cor <- names(which(abs(sapply(tumor_cors, function(x) x$vst[1,1])) < 0.7) )
# sapply(tumor_cors[bad_cor], function(x) x$vst_filt[1,1])

rda_pc10 <- lapply(tumor_pcs, function(y) rda(y$vst$x ~ y$scores$x))
# rda_pc2 <- lapply(tumor_pcs, function(y) rda(y$vst$x[, 1:2] ~ y$scores$x[, 1:2]))

var_prop10 <- sapply(rda_pc10, function(i) summary(i)$cont$importance[3, 10])
# var_prop2 <- sapply(rda_pc2, function(i) summary(i)$cont$importance[3, 2])

pc_prop10 <- sapply(tumor_pcs, function(x) cumsum(x$vst$sdev[1:10]**2)/sum(x$vst$sdev**2))

df_pcs_prad <- data.frame(PC1 = c(tumor_pcs[["TCGA-PRAD"]]$vst$x[, 1], tumor_pcs[["TCGA-PRAD"]]$scores$x[, 1]),
                     PC2 = c(tumor_pcs[["TCGA-PRAD"]]$vst$x[, 2], tumor_pcs[["TCGA-PRAD"]]$scores$x[, 2]),
                      dataset = rep(c("Original gene expression", "Gene set activities"), c(nrow(tumor_pcs[["TCGA-PRAD"]]$vst$x), nrow(tumor_pcs[["TCGA-PRAD"]]$scores$x))),
                    Sample = rep(factor(tumor_pcs[["TCGA-PRAD"]]$status), 2)) %>%
                    mutate(dataset = factor(dataset, levels = c("Original gene expression", "Gene set activities")),
                          Sample = ifelse(Sample == "Solid Tissue Normal", "Normal", "Cancer"),
                          Sample = factor(Sample, levels = c("Normal", "Cancer")))

pcs_model <- ggplot(df_pcs_prad, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_manual(values = c("green", "red")) +
ggtitle("PRAD dataset") +
  theme(plot.title = element_text(hjust = 0.5))

#
pcs_model_prad_ori <- data.frame(tumor_pcs[["TCGA-PRAD"]]$vst$x,
                       Sample = factor(tumor_pcs[["TCGA-PRAD"]]$status)) %>%
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
                       Sample = factor(tumor_pcs[["TCGA-PRAD"]]$status)) %>%
                    mutate( Sample = ifelse(Sample == "Solid Tissue Normal", "Normal", "Cancer"),
                            Sample = factor(Sample, levels = c("Normal", "Cancer"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Sample)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = c("green", "red")) +
    ggtitle("PRAD - NetActivity scores") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(sprintf("PC1 (%.1f%%)", tumor_pcs[["TCGA-PRAD"]]$scores$sdev[1]**2/sum(tumor_pcs[["TCGA-PRAD"]]$scores$sdev**2)*100)) +
    ylab(sprintf("PC2 (%.1f%%)", tumor_pcs[["TCGA-PRAD"]]$scores$sdev[2]**2/sum(tumor_pcs[["TCGA-PRAD"]]$scores$sdev**2)*100))



corplot_prad <- ggcorrplot(tumor_cors[["TCGA-PRAD"]], method = "circle", hc.order = FALSE,
      title = "PRAD") +
  scale_x_discrete("Original gene expression") +
  scale_y_discrete(name = "Gene set activities") +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
   axis.title.x = element_text(angle = 0, color = 'grey20', size = 20),
   axis.title.y = element_text(angle = 90, color = 'grey20', size = 20))



## IBD cohort
vst_ibd <- loadHDF5SummarizedExperiment("results/SRP042228/", prefix = "vsd_norm")

pc_vst_ibd <- prcomp(t(assay(vst_ibd)), rank. = 10)

assay(vst_ibd) <- data.matrix(assay(vst_ibd))
preproc_ibd <- prepareSummarizedExperiment(vst_ibd, "gtex_gokegg")
scores_ibd <- computeGeneSetScores(preproc_ibd, "gtex_gokegg")
pc_scores_ibd <- prcomp(t(assay(scores_ibd)), rank. = 10)
rda_ibd <- rda(pc_vst_ibd$x ~ pc_scores_ibd$x)

pc_summary <- data.frame(prop10 = c(var_prop10, summary(rda_ibd)$cont$importance[3, 10]),
          dataset = c(names(var_prop10), "IBD"),
          pc10 = c(pc_prop10[10, ], sum(pc_vst_ibd$sdev[1:10]**2)/sum(pc_vst_ibd$sdev**2)))

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
    xlab("Variance explained by the scores (%)") +
    ylab("Cumulative variance (%)") +
    ggtitle("PCs 1-10 correlation") +
      theme(plot.title = element_text(hjust = 0.5))

# png("figures/tumors_score_variance.png", width = 1300)
# plot_grid(pc_plot1, pc_plot10, ncol = 2)
# dev.off()

png("figures/tumors_score_variance.png")
pc_plot10
dev.off()


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


df_pcs_ibd <- data.frame(PC1 = c(pc_vst_ibd$x[, 1], -pc_scores_ibd$x[, 1]),
                     PC2 = c(pc_vst_ibd$x[, 2], -pc_scores_ibd$x[, 2]),
                      dataset = rep(c("Original gene expression", "Gene set activities"), c(nrow(pc_vst_ibd$x), nrow(pc_scores_ibd$x))),
                    Sample = rep(factor(vst_ibd$diagnosis), 2)) %>%
                    mutate(dataset = factor(dataset, levels = c("Original gene expression", "Gene set activities")),
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
pcs_model_ibd_scores <- data.frame(-pc_scores_ibd$x,
                       Sample = factor(vst_ibd$diagnosis)) %>%
                    mutate( Sample = recode(Sample, Control = "Normal"),
                            Sample = factor(Sample, levels = c("Normal", "CD", "UC"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Sample)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = c("green", "blue", "orange")) +
    ggtitle("IBD - NetActivity scores") +
    theme(plot.title = element_text(hjust = 0.5))  +
    xlab(sprintf("PC1 (%.1f%%)", pc_scores_ibd$sdev[1]**2/sum(pc_scores_ibd$sdev**2)*100)) +
    ylab(sprintf("PC2 (%.1f%%)", pc_scores_ibd$sdev[2]**2/sum(pc_scores_ibd$sdev**2)*100))



corplot_ibd <- ggcorrplot(cor(pc_vst_ibd$x, pc_scores_ibd$x), method = "circle", hc.order = FALSE,
      title = "IBD") +
  scale_x_discrete("Original gene expression") +
  scale_y_discrete(name = "Gene set activities") +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
   axis.title.x = element_text(angle = 0, color = 'grey20', size = 20),
   axis.title.y = element_text(angle = 90, color = 'grey20', size = 20))


#
png("figures/evaluationTCGAPRAD_panel2.png", width = 900, height = 1100)
plot_grid(plot_grid(pcs_model_prad_ori, pcs_model_prad_scores, rel_widths = c(6, 7)), plot_grid(pcs_model_ibd_ori, pcs_model_ibd_scores, rel_widths = c(6, 7)), plot_grid(corplot_prad, corplot_ibd, nrow = 1, labels = c("C", "D")),  ncol = 1, labels = c("A", "B", ""))
dev.off()
