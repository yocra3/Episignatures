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
--hdf5_dir data/GSE84727/ --input_probes results/preprocess/GEO_ref_blood/2021-06-11/cpg_medians.Rdata \
--trained_model results/GEOrefblood_model//2021-06-15/model.pb  -profile docker -resume

## Retrain model in GSE84727 data
nextflow run workflows/train_model_gse84727.nf  --hdf5_prefix GSE84727_raw \
--hdf5_dir data/GSE84727/ --network1 conf/network_params/network_2.py -profile docker -resume

## Apply model in GSE97362 data
nextflow run workflows/train_model_gse97362.nf  --hdf5_prefix GSE97362_raw \
--hdf5_dir data/GSE97362/ --input_probes results/preprocess/GEO_ref_blood/2021-06-11/cpg_medians.Rdata \
--trained_model results/GEOrefblood_model/2021-07-15/model.pb -profile docker -resume

nextflow run workflows/train_model_gse97362.nf  --hdf5_prefix GSE97362_raw \
--hdf5_dir data/GSE97362/ --input_probes --input_probes results/preprocess_tcga/2021-05-17/cpg_medians.Rdata \
--fname gse97362_tcga --trained_model results/tcga_model/2021-07-19/model.pb -profile docker -resume


## Retrain model in GSE97362 data
nextflow run workflows/train_model_gse97362.nf  --hdf5_prefix GSE97362_raw \
--hdf5_dir data/GSE97362/ --network1 conf/network_params/network_2.py -profile docker -resume

## Prepare GSE55763 files for training
nextflow run workflows/prepare_GSE55763_data.nf  --hdf5_prefix GSE55763_withNAs_noReplicates \
--hdf5_dir data/GSE55763/ -profile docker -resume

## Train GSE55763
nextflow run workflows/train_model_gse55763.nf  --hdf5_file results/preprocess/GSE55763/2021-06-08/assay_reshaped.h5 \
--config1 conf/network_params/random_search_1.py --config2 conf/network_params/random_search_2.py \
--config3 conf/network_params/random_search_3.py \
--network1 conf/network_params/network_3.py -profile docker -resume


## Prepare GEO ref blood files for training
nextflow run workflows/prepare_georefblood.nf  --hdf5_prefix GEOrefblood \
--hdf5_dir data/GEO_ref_blood/ -profile docker -resume

## Train GEO ref blood files
nextflow run workflows/train_model_georefblood.nf  --hdf5_file results/preprocess/GEO_ref_blood/2021-06-11/assay_reshaped.h5 \
--config1 conf/network_params/random_search_7.py \
--network1 conf/network_params/network_4.py -profile docker -resume
