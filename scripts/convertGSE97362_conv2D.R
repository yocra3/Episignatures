#'#################################################################################
#'#################################################################################
# Convert GSE97362 data to conv2D
#'#################################################################################
#'#################################################################################

## Load libraries
library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)
library(rhdf5)
library(limma)
library(pheatmap)
library(e1071)

## Load data
se <- loadHDF5SummarizedExperiment("data/GSE97362/", prefix = "GSE97362_raw")
se$disease <- factor(se$`disease state:ch1`, levels = c("Control", "CHARGE", "Kabuki",  "CHD7 variant","KMT2D variant"))
se$sex <- ifelse(se$characteristics_ch1 == "gender: female", "female", "male")
se$age <- as.numeric(gsub("age (years): ", "", se$characteristics_ch1.1, fixed = TRUE ))

## Sort SE
se <- se[!is.na(rowData(se)$MAPINFO)]
rowRanges(se) <- makeGRangesFromDataFrame(rowData(se), start.field = "MAPINFO", 
                                          end.field = "MAPINFO", ignore.strand = TRUE)
se <- sort(se)


## Prepare quantile
cont <- se[, se$disease == "Control"]
quant <- rowQuantiles(assay(cont), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)

imp_mat <- meffil:::impute.matrix(data.matrix(assay(se)), margin = 2)

imp_mat7 <- imp_mat[as.character(seqnames(se)) == 7, ]
new_array <- lapply(seq_len(ncol(imp_mat)), function(x) 
  cbind(imp_mat[, x], quant, start(se)))
new_array <- array(unlist(new_array), dim = c(nrow(imp_mat), ncol(quant) + 2, ncol(imp_mat7)))

discov_array <- new_array[, , ! se$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]
val_array <- new_array[, , se$`sample type:ch1` %in% c("Control for validation cohort",  "Validation cohort") ]

h5createFile("results/gse97362_episignature_validation/mini_conv_2D.h5")
h5write(discov_array, "results/gse97362_episignature_validation/mini_conv_2D.h5", "discov")
h5write(val_array, "results/gse97362_episignature_validation/mini_conv_2D.h5", "val")
h5closeAll()

# Check conv2D layers ####
discov <- se[, ! se$`sample type:ch1`  %in% c("Control for validation cohort",  "Validation cohort") ]
val <- se[,  se$`sample type:ch1` %in% c("Control for validation cohort",  "Validation cohort") ]
diseases <- c("CHARGE", "Kabuki")
names(diseases) <- c("CHARGE", "Kabuki")


layers <- "./results/gse97362_episignature_validation/conv_2D_outputs.h5" 
layers_train <- "./results/gse97362_episignature_validation/conv_2D_outputs_train.h5"

## Flatten ####
flatten <- h5read(layers, "flatten")

layers_train <- "./results/gse97362_episignature_validation/conv_2D_outputs_train.h5"
flatten_train <- h5read(layers_train, "flatten")

discov_flatten <- SummarizedExperiment(flatten_train, colData = colData(discov))
rownames(discov_flatten) <- paste0("Feat", seq_len(nrow(discov_flatten)))

val_flatten <- SummarizedExperiment(flatten, colData = colData(val))
rownames(val_flatten) <- paste0("Feat", seq_len(nrow(val_flatten)))

## Copy getFeaures from checkEpisignature.R

feats <- lapply(diseases, getFeatures, set = discov_flatten)



## Draw heatmaps ####
col_colors <- list(
  disease = c("Control" = "black", "CHARGE" = "green", "Kabuki" = "blue",
              "CHD7 variant" = "lightgreen", "KMT2D variant" = "cyan"),
  sex = c("female" = "purple", "male" = "lightblue")
)

## CHARGE
pheatmap(assay(discov_flatten)[feats$CHARGE[1:10000], discov_flatten$disease %in% c("Control", "CHARGE")], scale = "none", 
         annotation_col  = data.frame(colData(discov_flatten)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val_flatten)[feats$CHARGE[1:10000], val_flatten$disease %in% c("Control", "CHARGE", "CHD7 variant")], scale = "none", 
         annotation_col  = data.frame(colData(val_flatten)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)

## Kabuki
pheatmap(assay(discov_flatten)[feats$Kabuki[1:100], discov_flatten$disease %in% c("Control", "Kabuki")], scale = "none", 
         annotation_col  = data.frame(colData(discov_flatten)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val_flatten)[feats$Kabuki[1:100], val_flatten$disease %in% c("Control", "Kabuki", "KMT2D variant")], scale = "none", 
         annotation_col  = data.frame(colData(val_flatten)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)


## Dense1 ####
dense1 <- h5read(layers, "dense")
dense1_train <- h5read(layers_train, "dense")

discov_dense1 <- SummarizedExperiment(dense1_train, colData = colData(discov))
rownames(discov_dense1) <- paste0("Feat", seq_len(nrow(discov_dense1)))

val_dense1 <- SummarizedExperiment(dense1, colData = colData(val))
rownames(val_dense1) <- paste0("Feat", seq_len(nrow(val_dense1)))

feats_dense1 <- lapply(diseases, getFeatures, set = discov_dense1)

## Draw heatmaps ####
col_colors <- list(
  disease = c("Control" = "black", "CHARGE" = "green", "Kabuki" = "blue",
              "CHD7 variant" = "lightgreen", "KMT2D variant" = "cyan"),
  sex = c("female" = "purple", "male" = "lightblue")
)

## CHARGE
pheatmap(assay(discov_dense1)[feats_dense1$CHARGE, discov_dense1$disease %in% c("Control", "CHARGE")], scale = "none", 
         annotation_col  = data.frame(colData(discov_dense1)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val_dense1)[feats_dense1$CHARGE, val_dense1$disease %in% c("Control", "CHARGE", "CHD7 variant")], scale = "none", 
         annotation_col  = data.frame(colData(val_dense1)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
comb_dense1 <- cbind(val_dense1, discov_dense1)
pheatmap(assay(comb_dense1)[feats_dense1$CHARGE, comb_dense1$disease %in% c("Control", "CHARGE", "CHD7 variant")], scale = "none", 
         annotation_col  = data.frame(colData(comb_dense1)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)

## Kabuki
pheatmap(assay(discov_dense1)[feats_dense1$Kabuki, discov_dense1$disease %in% c("Control", "Kabuki")], scale = "none", 
         annotation_col  = data.frame(colData(discov_dense1)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(val_dense1)[feats_dense1$Kabuki, val_dense1$disease %in% c("Control", "Kabuki", "KMT2D variant")], scale = "none", 
         annotation_col  = data.frame(colData(val_dense1)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)
pheatmap(assay(comb_dense1)[feats_dense1$Kabuki, comb_dense1$disease %in% c("Control", "Kabuki", "KMT2D variant")], scale = "none", 
         annotation_col  = data.frame(colData(comb_dense1)[, c("disease", "sex")]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE)

## Copy trainSVm
svms <- lapply(diseases, trainSVM, set = discov_dense1, feats = feats_dense1)
table(predict(svms$CHARGE, t(assay(val_dense1))), val_dense1$disease)
table(predict(svms$Kabuki, t(assay(val_dense1))), val_dense1$disease)


## Explore kernels
val_array <- h5read("results/gse97362_episignature_validation/mini_conv_2D.h5", "val")
maxfeat <- which.max(rowSds(flatten))
cpg <- (maxfeat-1)* 4 + 1

df <- data.frame(start =  start(val)[cpg:(cpg+3)], assay(val[cpg:(cpg+3), ]))
df.cot <- gather(df, sample, methy, -1)
cot <- ifelse(flatten[39129, ] > 1, "Active", "Inactive")
names(cot) <- colnames(val)
df.cot$class <- cot[df.cot$sample]
ggplot(df.cot, aes(x = start, y = methy, col = class)) + geom_point()
ggplot(df.cot, aes(x = class, y = methy, col = class)) + geom_boxplot() +
  facet_grid(~ start)

sel <- val_array[cpg:(cpg+3), , ]
kernel <- matrix(c(0.01654436,  0.05999675,  0.04802385,  0.14613894,
                 -0.08279812, -0.16063476,  0.0498579 , -0.1944503 ,                                                                                                                                               
                  0.02151858,  0.00664616, -0.07187054,  0.06323695,
                 -0.03702331, -0.21502548,  0.02967073, -0.1081168 ,
                 -0.04659631,  0.19072166,  0.08757145,  0.0677878 ,
                  0.10209478,  0.18055128,  0.10955484, -0.18427055,
                 -0.19625418, -0.06793891, -0.00763006,  0.00837696), nrow = 7, 
                 byrow = TRUE)
vals <- sapply(seq_len(dim(sel)[3]) , function(x) {
  max(0, sum(kernel * t(sel[, , x])))
  })
  