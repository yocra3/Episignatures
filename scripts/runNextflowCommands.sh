#'#################################################################################
#'#################################################################################
# Nextflow commands
#'#################################################################################
#'#################################################################################
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiAzMzE2fS5kNjdiNWM3OGVkYWUyYWE3ZjM1OWYxOGNlYTU2NDBhYmMzNjdlODY2

## Prepare genotype files for imputation
nextflow run workflows/prepare_TCGA_data.nf  --hdf5_prefix methyTCGA \
--hdf5_dir data/tcga_hdf5/ -profile docker -resume
#nextflow run workflows/prepare_TCGA_data.nf  --hdf5_prefix data/tcga_hdf5/methyTCGA -profile docker -resume -with-tower
