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


## Generate final SE  ####
tcgaSE <-loadHDF5SummarizedExperiment(dir = "data/tcga_hdf5/tmp/", prefix = "tcga33_")
tcgaSE$project <- colData(tcgaSE) %>%
  data.frame() %>%
  left_join(mutate(query$results[[1]], barcode = cases) %>% 
              select(project, barcode), by = "barcode") %>%
  .$project
saveHDF5SummarizedExperiment(tcgaSE, "data/tcga_hdf5", prefix = "methyTCGA_")

## Filter CpGs with all missings  ####
pNA <- rowMeans(is.na(assay(tcgaSE)))
tcgaSE.filt <- tcgaSE[pNA < 0.9, ]
saveHDF5SummarizedExperiment(tcgaSE.filt, "data/tcga_hdf5", prefix = "methyTCGAfilt_")

## Create mini dataset for debugging ####
tcgaSE <-loadHDF5SummarizedExperiment(dir = "data/tcga_hdf5", prefix = "methyTCGAfilt_")
tcgaSE.mini <- tcgaSE[, c(1:10, 1001:1010)]
saveHDF5SummarizedExperiment(tcgaSE.mini, "data/tcga_hdf5", prefix = "testTCGA_")

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
  
  
# GSE90496 ####
## We do not download GSE109381 (superseries with also the validation cohort), because the methylation subtypes were
## only manually checked and defined in the reference dataset
library(GEOquery)
library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)

geo.full <- getGEO("GSE90496", GSEMatrix = TRUE)[[1]]
pheno <- pData(geo.full)[, c("methylation class:ch1", "characteristics_ch1", "geo_accession", "title")]
pheno$project <- pheno$`methylation class:ch1`
probes <- read_delim("data/GSE90496/probes.txt.gz", delim = "\t")

beta_raw <- read_delim("data/GSE90496/beta_vals.txt.gz", delim = "\t")
beta_raw <- data.matrix(beta_raw)
colnames(beta_raw) <- colnames(geo.full)
rownames(beta_raw) <- probes$ID_REF

gse90496 <-  SummarizedExperiment(DelayedArray(beta_raw), colData = pheno)
saveHDF5SummarizedExperiment(gse90496, "data/GSE90496/", prefix = "GSE90496_raw")


# GSE82273 ####
library(GEOquery)
library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)

options(timeout = 10000)

geo.full <- getGEO("GSE82273", GSEMatrix = TRUE)[[1]]
gse82273 <-  SummarizedExperiment(DelayedArray(exprs(geo.full)), 
                                  colData = pData(geo.full),
                                  rowData = fData(geo.full))
saveHDF5SummarizedExperiment(gse82273, "data/GSE82273/", prefix = "GSE82273_raw")

# GSE84727 ####
library(GEOquery)
library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)

options(timeout = 10000)

geo.full <- getGEO("GSE84727", GSEMatrix = TRUE)[[1]]

beta_raw <- read_delim("data/GSE84727/GSE84727_normalisedBetas.csv.gz", delim = ",")
beta_mat <- data.matrix(beta_raw[, -1])
rownames(beta_mat) <- beta_raw$X1


pheno <- pData(geo.full)
rownames(pheno) <- pheno$`sentrixids:ch1`
colnames(beta_mat) <- pheno[colnames(beta_mat), "geo_accession"]

pheno <- pData(geo.full)
pheno <- pheno[colnames(beta_mat), ]
pheno$project <- c("Control", "Schizophrenia")[as.numeric(pheno$`disease_status:ch1`)]

gse84727 <-  SummarizedExperiment(DelayedArray(beta_mat), 
                                  colData = pheno)
saveHDF5SummarizedExperiment(gse84727, "data/GSE84727/", prefix = "GSE84727_raw")

# GSE97362 ####
library(GEOquery)
library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)

options(timeout = 10000)

geo.full <- getGEO("GSE97362", GSEMatrix = TRUE)[[1]]

## Remove samples with VUS - some with control episignature, others disease episignature
geo.filt <- geo.full[, !grepl("VUS", geo.full$title)]
geo.filt$project <- geo.filt$characteristics_ch1.3GSE55763
geo.filt$project[geo.filt$project == "disease state: CHD7 variant"] <- "disease state: CHARGE"
geo.filt$project[geo.filt$project == "disease state: KMT2D variant"] <- "disease state: Kabuki"
geo.filt <- geo.filt[, geo.filt$characteristics_ch1.3 != "disease state: KDM6A variant"]
gse97362 <-  SummarizedExperiment(DelayedArray(exprs(geo.filt)), 
                                  colData = pData(geo.filt),
                                  rowData = fData(geo.filt))
saveHDF5SummarizedExperiment(gse97362, "data/GSE97362/", prefix = "GSE97362_raw")

# GSE55763 ####
library(GEOquery)
library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


options(timeout = 10000)

geo.full <- getGEO("GSE55763", GSEMatrix = TRUE)[[1]]

pheno <- pData(geo.full)[, c("characteristics_ch1.1", "geo_accession", "title", 
                             "characteristics_ch1.2", "characteristics_ch1.3", 
                             "age:ch1", "dataset:ch1",  "gender:ch1")]
pheno$project <- pheno$`gender:ch1`
pheno_cols <- gsub("Peripheral blood, ", "", pheno$title)

probes <- read_delim("data/GSE55763/probes.txt.gz", delim = "\t")

beta_raw <- read_delim("data/GSE55763/beta_vals.txt.gz", delim = "\t")
beta_raw <- data.matrix(beta_raw[, pheno_cols])

colnames(beta_raw) <- colnames(geo.full)
rownames(beta_raw) <- probes$ID_REF

gse55763 <-  makeG(DelayedArray(beta_raw), colData = pheno)
saveHDF5SummarizedExperiment(gse55763, "data/GSE55763/", prefix = "GSE55763_raw")

## Add NA to positions with bad detection p-values
se <- loadHDF5SummarizedExperiment(dir = "data/GSE55763/", prefix = "GSE55763_raw")
cols <- read_delim("data/GSE55763/betas_header.txt", delim = "\t")
detP <- read_delim("data/GSE55763/detection_pvals.txt.gz", delim = "\t")
colnames(detP) <- colnames(cols)
arr_names <- gsub("Peripheral blood, ", "", se$title)
detP <- DelayedArray(data.matrix(detP[, arr_names]))
assay(se)[detP != 0] <- NA

## Add annotation
data("Locations")
grCpG <- makeGRangesFromDataFrame(Locations, start.field = "pos", end.field = "pos")
rowRanges(se) <- grCpG[rownames(se)]
saveHDF5SummarizedExperiment(se, "data/GSE55763/", prefix = "GSE55763_withNAs")

se <- loadHDF5SummarizedExperiment(dir = "data/GSE55763/", prefix = "GSE55763_withNAs")
se <- se[, se$characteristics_ch1.1 == "dataset: population study"]
saveHDF5SummarizedExperiment(se, dir = "data/GSE55763/", prefix = "GSE55763_withNAs_noReplicates_")

# Create GEO blood reference ####
library(HDF5Array)
library(SummarizedExperiment)

se1 <- loadHDF5SummarizedExperiment("data/GSE82273/", prefix = "GSE82273_raw")
se2 <- loadHDF5SummarizedExperiment(dir = "data/GSE55763/", prefix = "GSE55763_withNAs_noReplicates_")

se1$sex <- se1$`gender:ch1`
se2$sex <- ifelse(se2$project == "F", "female", "male")

colData(se1) <- colData(se1)[, c("title", "geo_accession", "sex")]
colData(se2) <- colData(se2)[, c("title", "geo_accession", "sex")]

se_comb <- cbind(se1[rownames(se2), ], se2)
se_comb$project <- se_comb$sex
rowRanges(se_comb) <- rowRanges(se2)
saveHDF5SummarizedExperiment(se_comb, dir = "data/GEO_ref_blood/", prefix = "GEOrefblood_")
