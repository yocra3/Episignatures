#'#################################################################################
#'#################################################################################
# Nextflow commands
#'#################################################################################
#'#################################################################################
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiAzMzE2fS5kNjdiNWM3OGVkYWUyYWE3ZjM1OWYxOGNlYTU2NDBhYmMzNjdlODY2

## Prepare TCGA files for imputation
nextflow run workflows/prepare_TCGA_data.nf  --hdf5_prefix methyTCGA \
--hdf5_dir data/tcga_hdf5/ -profile docker -resume

## Train model in TCGA data
nextflow run workflows/train_model_TCGA.nf  --hdf5_file results/preprocess/2021-05-03/assay_reshaped.h5 \
--config1 conf/network_params/random_search_1.py --network1 conf/network_params/network_1.py -profile docker -resume
