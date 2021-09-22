#'#################################################################################
#'#################################################################################
# Check GSE97362 episignature
#'#################################################################################
#'#################################################################################

library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(pheatmap)
library(e1071)

se <- loadHDF5SummarizedExperiment("data/GSE97362/", prefix = "GSE97362_raw")
se$disease <- factor(se$`disease state:ch1`, levels = c("Control", "CHARGE", "Kabuki",  "CHD7 variant","KMT2D variant"))
se$sex <- ifelse(se$characteristics_ch1 == "gender: female", "female", "male")
se$age <- as.numeric(gsub("age (years): ", "", se$characteristics_ch1.1, fixed = TRUE ))

discov <- se[, ! se$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]
val <- se[, se$`sample type:ch1` %in% c("Control for validation cohort",  "Validation cohort") ]



diseases <- c("CHARGE", "Kabuki")
names(diseases) <- c("CHARGE", "Kabuki")

getFeatures <- function(set, disease){
  
  set <- set[, set$disease %in% c("Control", disease)]
  set$disease <- droplevels(set$disease)
  model <- model.matrix(~ disease + sex + age, colData(set))
  
  lmF <- lmFit(assay(set[, rownames(model)]), model)
  lmFe <- eBayes(lmF)
  tabA <- topTable(lmFe, coef = 2, n = Inf)
  selTab <- subset(tabA, !is.na(logFC) & adj.P.Val < 0.01 & abs(logFC) > 0.1)
  selCpGs <- rownames(selTab)
  
}
feats <- lapply(diseases, getFeatures, set = discov)

featsDF <- data.frame(cpg = unlist(feats), disease = rep(names(feats), lengths(feats)))
write.table(featsDF, file = "results/gse97362_episignature_validation/svm_paper_feats.txt", quote = FALSE,
            row.names = FALSE)


discov_all <- discov
assay(discov_all) <- meffil:::impute.matrix(data.matrix(assay(discov_all)), margin = 2)
saveHDF5SummarizedExperiment(discov_all, "results/gse97362_episignature_validation/", prefix = "discov_")

discov_feats <- discov[unique(unlist(feats)), ]
assay(discov_feats) <- meffil:::impute.matrix(data.matrix(assay(discov_feats)), margin = 2)

saveHDF5SummarizedExperiment(discov_feats, "results/gse97362_episignature_validation/", prefix = "discov_allFeats")
saveHDF5SummarizedExperiment(discov_feats[feats$CHARGE, ], "results/gse97362_episignature_validation/", prefix = "discov_charge")
saveHDF5SummarizedExperiment(discov_feats[feats$Kabuki, ], "results/gse97362_episignature_validation/", prefix = "discov_kabuki")
write.table(discov$disease, file = "results/gse97362_episignature_validation/discov_labels.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

val_all <- val
assay(val_all) <- meffil:::impute.matrix(data.matrix(assay(val_all)), margin = 2)
saveHDF5SummarizedExperiment(val_all, "results/gse97362_episignature_validation/", prefix = "val_")



val_feats <- val[unique(unlist(feats)), ]
assay(val_feats) <- meffil:::impute.matrix(data.matrix(assay(val_feats)), margin = 2)

saveHDF5SummarizedExperiment(val_feats, "results/gse97362_episignature_validation/", prefix = "val_allFeats")
saveHDF5SummarizedExperiment(val_feats[feats$CHARGE, ], "results/gse97362_episignature_validation/", prefix = "val_charge")
saveHDF5SummarizedExperiment(val_feats[feats$Kabuki, ], "results/gse97362_episignature_validation/", prefix = "val_kabuki")
write.table(val$disease, file = "results/gse97362_episignature_validation/val_labels.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

## Draw heatmaps ####
col_colors <- list(
  disease = c("Control" = "black", "CHARGE" = "green", "Kabuki" = "blue",
              "CHD7 variant" = "lightgreen", "KMT2D variant" = "cyan"),
  sex = c("female" = "purple", "male" = "lightblue")
)

## CHARGE
pheatmap(assay(discov)[feats$CHARGE, discov$disease %in% c("Control", "CHARGE")], scale = "none", 
         annotation_col  = data.frame(colData(discov)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val)[feats$CHARGE, val$disease %in% c("Control", "CHARGE", "CHD7 variant")], scale = "none", 
         annotation_col  = data.frame(colData(val)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)

## Kabuki
pheatmap(assay(discov)[feats$Kabuki, discov$disease %in% c("Control", "Kabuki")], scale = "none", 
         annotation_col  = data.frame(colData(discov)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val)[feats$Kabuki, val$disease %in% c("Control", "Kabuki", "KMT2D variant")], scale = "none", 
         annotation_col  = data.frame(colData(val)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)

## Train SVM ####
trainSVM <-  function(set, disease, feats){
  set <- set[, set$disease %in% c("Control", disease)]
  selCpGs <- feats[[disease]]
  mat <- data.matrix(assay(set[selCpGs, ]))
  imp_mat <- meffil:::impute.matrix(mat, margin = 1, fun = function(x) median(x, na.rm = T))
  
  df <- data.frame(pathClass = factor(set$disease), t(imp_mat))
  model_svm <- svm(pathClass ~ ., df)
}
svms <- lapply(diseases, trainSVM, set = discov, feats = feats)

all_mat <- data.matrix(t(assay(se[unique(unlist(feats)), ])))
imp_mat <- meffil:::impute.matrix(all_mat, margin = 2, fun = function(x) median(x, na.rm = T))

pred_charge <- predict(svms$CHARGE, imp_mat)
table(pred_charge, se$`disease state:ch1`)

pred_kabuki <- predict(svms$Kabuki, imp_mat)
table(pred_kabuki, se$`disease state:ch1`)

## Overlap con referencia ####
load("results/preprocess/GEO_ref_blood/2021-06-11/cpg_medians.Rdata")
mean(unique(featsDF$cpg) %in% names(medians))

medians_feats <- rowMedians(data.matrix(assay(se[unique(featsDF$cpg), ])), na.rm = TRUE)
names(medians_feats) <- unique(featsDF$cpg)

com <- intersect(names(medians_feats) , names(medians))
plot(medians_feats[com], medians[com])

## Testear en datos modificados ####
library(rhdf5)
file <- "./work/50/b41517a3a23b61c1061f8991c67e1b/assay_reshaped.h5" ## Link to nextflow folder
h5mat <- h5read(file, "methy")[1, , ]
rownames(h5mat) <- names(medians)
colnames(h5mat) <- colnames(se)
h5mat <- h5mat[com, ]
out_probes <- setdiff(names(medians_feats), rownames(h5mat))
out <- matrix(rep(medians_feats[out_probes], ncol(h5mat)), 
                           nrow = length(out_probes), ncol = ncol(h5mat),
                           dimnames = list(out_probes, colnames(h5mat)))
h5mat_f <- rbind(h5mat, out)
table(predict(svms$CHARGE, t(h5mat_f)), se$`disease state:ch1`)

table(predict(svms$Kabuki, t(h5mat_f)), se$`disease state:ch1`)
table(pred_kabuki, se$`disease state:ch1`)


# Comprobar features ####
## Flatten layer ####
flatten <- read_delim("results/gse97362_model/2021-07-16/flatten_1.tsv", delim = "\t")
flatten <- data.matrix(flatten)

se_flatten <- SummarizedExperiment(t(flatten), colData = colData(se))
rownames(se_flatten) <- paste0("Features", seq_len(nrow(se_flatten)))

discov_flatten <- se_flatten[, ! se_flatten$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]
val_flatten <- se_flatten[, se_flatten$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]

flatten_feats <- lapply(diseases, getFeatures, set = discov_flatten)


## CHARGE
pheatmap(assay(discov_flatten)[flatten_feats$CHARGE, discov_flatten$disease %in% c("Control", "CHARGE")], scale = "none", 
         annotation_col  = data.frame(colData(discov_flatten)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val_flatten)[flatten_feats$CHARGE, val_flatten$disease %in% c("Control", "CHARGE", "CHD7 variant")], scale = "none", 
         annotation_col  = data.frame(colData(val_flatten)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)

## Kabuki
pheatmap(assay(discov_flatten)[flatten_feats$Kabuki, discov_flatten$disease %in% c("Control", "Kabuki")], scale = "none", 
         annotation_col  = data.frame(colData(discov_flatten)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val_flatten)[flatten_feats$Kabuki, val_flatten$disease %in% c("Control", "Kabuki", "KMT2D variant")], scale = "none", 
         annotation_col  = data.frame(colData(val_flatten)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)


svms_flatten <- lapply(diseases, trainSVM, set = discov_flatten, feats = flatten_feats)

table(predict(svms_flatten$CHARGE, t(assay(se_flatten))), se_flatten$disease)
table(predict(svms_flatten$Kabuki, t(assay(se_flatten))), se_flatten$disease)

saveHDF5SummarizedExperiment(discov_flatten, "results/gse97362_episignature_validation/", prefix = "discov_wholeFlatten")
saveHDF5SummarizedExperiment(discov_flatten[unique(unlist(flatten_feats)), ], "results/gse97362_episignature_validation/", prefix = "discov_allFlatten")
saveHDF5SummarizedExperiment(discov_flatten[flatten_feats$CHARGE, ], "results/gse97362_episignature_validation/", prefix = "discov_chargeFlatten")
saveHDF5SummarizedExperiment(discov_flatten[flatten_feats$Kabuki, ], "results/gse97362_episignature_validation/", prefix = "discov_kabukiFlatten")

saveHDF5SummarizedExperiment(val_flatten, "results/gse97362_episignature_validation/", prefix = "val_wholeFlatten")
saveHDF5SummarizedExperiment(val_flatten[unique(unlist(flatten_feats)), ], "results/gse97362_episignature_validation/", prefix = "val_allFlatten")
saveHDF5SummarizedExperiment(val_flatten[flatten_feats$CHARGE, ], "results/gse97362_episignature_validation/", prefix = "val_chargeFlatten")
saveHDF5SummarizedExperiment(val_flatten[flatten_feats$Kabuki, ], "results/gse97362_episignature_validation/", prefix = "val_kabukiFlatten")


## Dense layer ####
dense <- read_delim("results/gse97362_model/2021-07-16/dense_1.tsv", delim = "\t")
dense <- data.matrix(dense)

se_dense <- SummarizedExperiment(t(dense), colData = colData(se))
rownames(se_dense) <- paste0("Dense", seq_len(nrow(se_dense)))

discov_dense <- se_dense[, ! se_dense$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]
val_dense <- se_dense[, se_dense$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]

dense_feats <- lapply(diseases, getFeatures, set = discov_dense)

## Kabuki
pheatmap(assay(discov_dense)[dense_feats$Kabuki, discov_dense$disease %in% c("Control", "Kabuki")], scale = "none", 
         annotation_col  = data.frame(colData(discov_dense)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val_dense)[dense_feats$Kabuki, val_dense$disease %in% c("Control", "Kabuki", "KMT2D variant")], scale = "none", 
         annotation_col  = data.frame(colData(val_dense)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
svm_dense <- trainSVM("Kabuki", set = discov_dense, feats = dense_feats)
table(predict(svm_dense, t(assay(se_dense))), se_dense$disease)


pheatmap(assay(se_dense), scale = "none", 
         annotation_col  = data.frame(colData(se_dense)[, c("disease", "sex"), drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)


## Flatten layer TCGA ####
flatten_tcga <- read_delim("results/gse97362_tcga/2021-07-19/flatten_1.tsv", delim = "\t")
flatten_tcga <- data.matrix(flatten_tcga)

se_flatten_tcga <- SummarizedExperiment(t(flatten_tcga), colData = colData(se))
rownames(se_flatten_tcga) <- paste0("Features", seq_len(nrow(se_flatten_tcga)))

discov_flatten_tcga <- se_flatten_tcga[, ! se_flatten_tcga$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]
val_flatten_tcga <- se_flatten_tcga[, se_flatten_tcga$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]

flatten_feats_tcga <- lapply(diseases, getFeatures, set = discov_flatten_tcga)


## Kabuki
pheatmap(assay(discov_flatten_tcga)[flatten_feats_tcga$Kabuki, discov_flatten_tcga$disease %in% c("Control", "Kabuki")], scale = "none", 
         annotation_col  = data.frame(colData(discov_flatten_tcga)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val_flatten_tcga)[flatten_feats_tcga$Kabuki, val_flatten_tcga$disease %in% c("Control", "Kabuki", "KMT2D variant")], scale = "none", 
         annotation_col  = data.frame(colData(val_flatten_tcga)[, "disease", drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)


svm_flatten_tcga <- trainSVM("Kabuki", set = discov_flatten_tcga, feats = flatten_feats_tcga)
table(predict(svm_flatten_tcga, t(assay(se_flatten_tcga))), se_flatten_tcga$disease)

## Dense layer tcga ####
dense_tcga <- read_delim("results/gse97362_tcga/2021-07-19/dense_1.tsv", delim = "\t")
dense_tcga <- data.matrix(dense_tcga)

se_dense_tcga <- SummarizedExperiment(t(dense_tcga), colData = colData(se))
rownames(se_dense_tcga) <- paste0("Dense", seq_len(nrow(se_dense_tcga)))

discov_dense_tcga <- se_dense_tcga[, ! se_dense_tcga$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]
val_dense_tcga <- se_dense_tcga[, se_dense_tcga$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]

dense_feats_tcga <- lapply(diseases, getFeatures, set = discov_dense_tcga)
pheatmap(assay(se_dense_tcga)[rowSums(assay(se_dense_tcga)) >0, ], scale = "none", 
         annotation_col  = data.frame(colData(se_dense_tcga)[, c("disease", "sex"), drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)

se_dense_tcga$comb <- ifelse(se_dense_tcga$disease == "CHD7 variant", "CHARGE", 
                             ifelse(se_dense_tcga$disease == "KMT2D variant", "Kabuki", 
                             as.character(se_dense_tcga$disease)))
mod <- model.matrix(~ comb - 1, colData(se_dense_tcga))
cors <- cor(t(assay(se_dense_tcga)), mod)
cors <- cors[rowMeans(is.na(cors)) != 1, ]

selcors <- cors[apply(cors, 1, function(x) any(abs(x) > 0.2)), ]

pheatmap(assay(se_dense_tcga)[rownames(selcors), ], scale = "none", 
         annotation_col  = data.frame(colData(se_dense_tcga)[, c("disease", "sex"), drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)

df <- data.frame(pathClass = factor(discov_dense_tcga$disease), t(assay(discov_dense_tcga)[rowSums(assay(discov_dense_tcga)) >0, ]))
model_svm <- svm(pathClass ~ ., df)
table(predict(model_svm, t(assay(se_dense_tcga))), se_dense_tcga$disease)


