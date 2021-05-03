#'#################################################################################
#'#################################################################################
# Download data
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(GEOquery)
library(TCGAbiolinks)
library(HDF5Array)
library(SummarizedExperiment)

# TCGA ####
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grep("TCGA", projects)]
query <- GDCquery(project = projects,
                  data.category = "DNA methylation",
                  legacy = TRUE,
                  platform = c("Illumina Human Methylation 450")
)
GDCdownload(query, method = "api", files.per.chunk = 10)

queryList <- lapply(projects, GDCquery, data.category = "DNA methylation",
                  legacy = TRUE, platform = c("Illumina Human Methylation 450"))
names(queryList) <- projects

selVars <- c("barcode", "patient", "sample", "shortLetterCode", "definition", 
             "sample_type_id", "sample_type", "tissue_type",
             "tumor_stage", "tissue_or_organ_of_origin", "days_to_last_follow_up", 
             "primary_diagnosis", "age_at_diagnosis", "prior_malignancy", 
             "classification_of_tumor", "tumor_grade", "gender", "ethnicity", 
             "vital_status", "race", "age_at_index", "days_to_death", "primary_site")

makeHDF5 <- function(query, name){
  message(name)
  if (file.exists(paste0("data/tcga_hdf5/", name, "_se.rds"))){
    return(NULL)
  }
  obj <- GDCprepare(query)
  assay(obj) <-  DelayedArray(assay(obj))
  colData(obj) <- colData(obj)[, selVars]
  saveHDF5SummarizedExperiment(obj, "data/tcga_hdf5", prefix= paste0(name, "_"))
  NULL
}
tcgaL <- lapply(projects, function(x) makeHDF5(queryList[[x]], x))

## Merge hdf5
prefixes <- gsub("assays.h5", "", dir("data/tcga_hdf5/", pattern = "h5"))
a <- lapply(prefixes, function(x){
  i <- which(prefixes == x)
  if (file.exists(paste0("data/tcga_hdf5/tmp/tcga", i, "_se.rds"))){
    return(NULL)
  }
  se <- loadHDF5SummarizedExperiment(dir = "data/tcga_hdf5/", prefix = x)
  if(i > 1){
    old <- loadHDF5SummarizedExperiment(dir = "data/tcga_hdf5/tmp", 
                                        prefix = paste0("tcga", i - 1, "_"))
    se <- cbind(se, old)
  }
  saveHDF5SummarizedExperiment(se, "data/tcga_hdf5/tmp", 
                               prefix = paste0("tcga", i, "_"))
} )


## Generate final SE
tcgaSE <-loadHDF5SummarizedExperiment(dir = "data/tcga_hdf5/tmp/", prefix = "tcga33_")
tcgaSE$project <- colData(tcgaSE) %>%
  data.frame() %>%
  left_join(mutate(query$results[[1]], barcode = cases) %>% 
              select(project, barcode), by = "barcode") %>%
  .$project
saveHDF5SummarizedExperiment(tcgaSE, "data/tcga_hdf5", prefix = "methyTCGA_")

## Filter CpGs with all missings
pNA <- rowMeans(is.na(assay(tcgaSE)))
tcgaSE.filt <- tcgaSE[pNA < 0.9, ]
saveHDF5SummarizedExperiment(tcgaSE.filt, "data/tcga_hdf5", prefix = "methyTCGAfilt_")


## Play with data
library(DelayedMatrixStats)
library(BiocSingular)
library(BiocParallel)
library(tidyverse) 
library(HDF5Array)
library(SummarizedExperiment)

tcga <- loadHDF5SummarizedExperiment("data/tcga_hdf5", prefix = "methyTCGAfilt_")
proj <- tcga$project
proj[tcga$sample_type == "Solid Tissue Normal"] <- "Normal"
write.table(proj, file = "data/tcga_hdf5/methyTCGAfilt_proj.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

tcga.full <- tcga

tcga.full <- tcga[pNA[pNA < 0.9] == 0, ]
pcs.out <- runPCA(t(assay(tcga.full)), center = FALSE, rank = 10, BSPARAM = FastAutoParam(),
                  BPPARAM = MulticoreParam(10))
pc.comb <- pcs.out$x %>%
  data.frame() %>%
  mutate(barcode = rownames(.)) %>%
  as_tibble() %>%
  left_join(colData(tcga.full) %>% data.frame(), by = "barcode") %>%
  left_join(mutate(query$results[[1]], sample = sample.submitter_id) %>% 
              select(project, sample), by = "sample") 

ggplot(pc.comb, aes(x = PC1, y = PC2, color = project)) +
  geom_point() +
  theme_bw()

ggplot(pc.comb, aes(x = PC1, color = sample_type)) +
  geom_density() +
  theme_bw()
  
  
# GSE109381 ####
pheno <- getGEO("GSE109381", GSEMatrix = TRUE)
# files <- getGEOSuppFiles("GSE109381", baseDir = "data/") ## Download data (manually download and decompress)

gzfiles <- dir("data/GEO/GSE62629", pattern = "gz", full.names = TRUE)

methy <- lapply(gzfiles, function(x) read.delim(gzfile(x), header = FALSE, as.is = TRUE))
methmat <- Reduce(cbind, lapply(methy, `[`, 3))
sampNames <- gsub("data/GEO/GSE62629/", "", gzfiles)
colnames(methmat) <- gsub("_.*", "", sampNames)

annot <- methy[[1]][, 1:2]
annot$start <- annot$V2 + 1 ## Add 1 to ensure matching with EPIC annotation
annotGR <- makeGRangesFromDataFrame(annot, seqnames.field = "V1",
                                    end.field = "start")
