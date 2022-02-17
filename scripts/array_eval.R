## Explore array
library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)

se <- loadHDF5SummarizedExperiment("data/GSE97362/", prefix = "GSE97362_raw")
se$disease <- factor(se$`disease state:ch1`, levels = c("Control", "CHARGE", "Kabuki",  "CHD7 variant","KMT2D variant"))
se$sex <- ifelse(se$characteristics_ch1 == "gender: female", "female", "male")
se$age <- as.numeric(gsub("age (years): ", "", se$characteristics_ch1.1, fixed = TRUE ))


se <- se[rowData(se)$CHR != "", ]
methy <- se[grep("cg", rownames(se)), ]
ranges <- GRanges(paste0(rowData(methy)$CHR, ":", rowData(methy)$AddressA_ID))
names(ranges) <- rownames(methy)
rowRanges(methy) <- ranges

methy <- sort(methy)

dist <- start(methy)[-1] -  start(methy)[-nrow(methy)] 
dist[dist < 0] <- 1e7

chr22 <- methy[seqnames(methy) == 22, ]
cor22 <- cor(t(data.matrix(assay(chr22))), use = "pairwise.complete.obs")
cor22_df <- data.frame(cor22)
cor22_df$CpG1 <- rownames(cor22_df)
cor22_tib <- gather(cor22_df, CpG2, cor, -ncol(cor22_df))
cor22_tib$CordDist <- start(rowRanges(chr22)[cor22_tib$CpG2]) - start(rowRanges(chr22)[cor22_tib$CpG1]) 

posdf <- data.frame(pos = seq_len(nrow(chr22)))
rownames(posdf) <- rownames(chr22)
cor22_tib$PosDist <- posdf[cor22_tib$CpG2, "pos"] - posdf[cor22_tib$CpG1, "pos"]
cor22_tib_filt <- filter(cor22_tib, CordDist > 0)
# ggplot(cor22_tib_filt, aes(x = CordDist, y = cor)) + geom_point() + theme_bw()
# cor22_tib_filt %>%
#   filter(CordDist < 4e5) %>%
# ggplot(aes(x = CordDist, y = cor)) + geom_point(alpha = 0.01) + theme_bw() 
cor22_tib_filt %>%
  mutate(DistBin = cut(CordDist, breaks = c(0, 1e3, 3e3, 10e3, 30e3, 100e3, 500e3, 1e9),
                       labels = c("<1Kb", "1-3 Kb", "3-10Kb", "10-30Kb", "30-100Kb", "100-500Kb", ">500Kb"))) %>%
  ggplot(aes(x = cor, color = DistBin)) + geom_density() + theme_bw() + geom_vline(xintercept = c(-0.3, 0, 0.3))
cor22_tib_filt %>%
  ggplot(aes(x = cor)) + geom_density() + theme_bw() 


cor22_tib_filt <- cor22_tib_filt %>%
  mutate(DistBin = cut(CordDist, breaks = c(0, 1e3, 3e3, 10e3, 30e3, 100e3, 500e3, 1e9),
                       labels = c("<1Kb", "1-3 Kb", "3-10Kb", "10-30Kb", "30-100Kb", "100-500Kb", ">500Kb")),
         PosBin = cut(PosDist, breaks = c(0, 5, 10, 20, 40, 60, 80, 300),
                      labels = c("1-5", "6-10", "11-20", "21-40", "41-60", "61-80", ">80")))
cor22_tib_filt %>%
  ggplot(aes(y = abs(cor), x = DistBin)) + geom_boxplot() + theme_bw()

prop.table(table(cor22_tib_filt$DistBin, abs(cor22_tib_filt$cor) > 0.3), margin = 1)

cor22_tib_filt %>%
  ggplot(aes(x = cor, color = PosBin)) + geom_density() + theme_bw() 

cor22_tib_filt %>%
  ggplot(aes(y = abs(cor), x = PosBin)) + geom_boxplot() + theme_bw()
prop.table(table(cor22_tib_filt$PosBin, abs(cor22_tib_filt$cor) > 0.3), margin = 1)


range22 <- apply(data.matrix(assay(chr22)), 1, function(x) quantile(x, 0.99, na.rm = TRUE) - quantile(x, 0.01, na.rm = TRUE))
cor22_tib_filt$range1 <- range22[cor22_tib_filt$CpG1]
cor22_tib_filt$range2 <- range22[cor22_tib_filt$CpG2]

cor22_tib_filt %>%
  mutate(range1_cat = cut(range1, breaks = c(0, 0.1, 0.2, 1),
                       labels = c("<0.1", "0.1-0.2", ">0.2")),
         range2_cat = cut(range2, breaks = c(0, 0.1, 0.2, 1),
                          labels = c("<0.1", "0.1-0.2", ">0.2"))) %>%
  ggplot(aes(x = cor, color = DistBin)) + geom_density() + theme_bw() + geom_vline(xintercept = c(-0.3, 0, 0.3)) +
  facet_grid(range1_cat ~ range2_cat)

cor22_tib_var <- cor22_tib_filt %>% 
  filter(range1 > 0.2 & range2 > 0.2)
prop.table(table(cor22_tib_var$DistBin, abs(cor22_tib_var$cor) > 0.3), margin = 1)
prop.table(table(cor22_tib_filt$DistBin, abs(cor22_tib_filt$cor) > 0.3, 
           ifelse(cor22_tib_filt$range1 < 0.2, 
                  ifelse(cor22_tib_filt$range2 < 0.2, "Both no Var", "No var1"), 
                  ifelse(cor22_tib_filt$range2 < 0.2, "No var2", "Both var"))), margin = 1)
