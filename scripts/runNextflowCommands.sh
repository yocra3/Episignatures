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
nextflow run workflows/train_model_TCGA.nf  --hdf5_file results/preprocess/2021-05-17/assay_reshaped.h5 \
--config1 conf/network_params/random_search_1.py --config2 conf/network_params/random_search_2.py \
--config3 conf/network_params/random_search_3.py \
--network1 conf/network_params/network_2.py -profile docker -resume

## Apply model in GSE90498 data
nextflow run workflows/train_model_gse90496.nf  --hdf5_prefix GSE90496_raw \
--hdf5_dir data/GSE90496/ --input_probes results/preprocess_tcga/2021-05-17/cpg_medians.Rdata \
--trained_model results/tcga_model/2021-05-27/model.pb -profile docker -resume

## Apply model in GSE84727 data
nextflow run workflows/train_model_gse84727.nf  --hdf5_prefix GSE84727_raw \
--hdf5_dir data/GSE84727/ --input_probes results/preprocess_tcga/2021-05-17/cpg_medians.Rdata \
--trained_model results/tcga_model/2021-05-27/model.pb -profile docker -resume

## Retrain model in GSE84727 data
nextflow run workflows/train_model_gse84727.nf  --hdf5_prefix GSE84727_raw \
--hdf5_dir data/GSE84727/ --network1 conf/network_params/network_2.py -profile docker -resume

## Apply model in GSE97362 data
nextflow run workflows/train_model_gse97362.nf  --hdf5_prefix GSE97362_raw \
--hdf5_dir data/GSE97362/ --input_probes results/preprocess_tcga/2021-05-17/cpg_medians.Rdata \
--trained_model results/tcga_model/2021-05-27/model.pb -profile docker -resume

nextflow run workflows/train_model_gse97362.nf  --hdf5_prefix GSE97362_filter \
--hdf5_dir data/GSE97362/ --input_probes results/preprocess_tcga/2021-05-17/cpg_medians.Rdata \
--trained_model results/tcga_model/2021-05-27/model.pb -profile docker -resume


## Retrain model in GSE97362 data
nextflow run workflows/train_model_gse97362.nf  --hdf5_prefix GSE97362_raw \
--hdf5_dir data/GSE97362/ --network1 conf/network_params/network_2.py -profile docker -resume
