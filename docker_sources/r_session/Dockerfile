# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_15

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("Gviz", "a4Base", "tidyverse", "limma", "e1071", "sva",\
  "cowplot", "dendextend", "DelayedMatrixStats", "BiocSingular", "BiocParallel", "HDF5Array",\
  "SummarizedExperiment", "hipathia", "fgsea", "pheatmap", "DESeq2", "reactome.db","GSVA"))'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
