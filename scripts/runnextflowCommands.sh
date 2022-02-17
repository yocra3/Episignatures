#'#################################################################################
#'#################################################################################
# Nextflow commands
#'#################################################################################
#'#################################################################################
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiAzMzE2fS5kNjdiNWM3OGVkYWUyYWE3ZjM1OWYxOGNlYTU2NDBhYmMzNjdlODY2

## Train model for blood
nextflow run workflows/train_model.nf --hdf5_file results/preprocess/GEO_ref_blood/2021-06-11/assay_reshaped.h5 \
--name blood_DNN1 --step random --network conf/network_params/dnn_network1.py \
--config conf/network_params/dnn1_random1.py -profile docker -resume

## Get features results for CNN models from TCGA and blood data
nextflow run workflows/train_model.nf --hdf5_file results/preprocess/gse97362/2021-06-02/assay_reshaped.h5 \
--name gse97362_TCGA --step features --model results/tcga_model/2021-05-27/model.pb -profile docker -resume

nextflow run workflows/train_model.nf --hdf5_file results/preprocess/gse97362/2021-07-16/assay_reshaped.h5 \
--name gse97362_Blood --step features --model  results/GEOrefblood_model/2021-07-15/model.pb -profile docker -resume

## Transfer learning
nextflow run workflows/train_model.nf --hdf5_file results/preprocess/gse97362/2021-06-02/assay_reshaped.h5 \
--name gse97362_TCGA --step transfer --layer "flatten_1" --model results/tcga_model/2021-05-27/model.pb -profile docker -resume


## Add distance to next CpG to methylation matrix
nextflow run workflows/train_model.nf --hdf5_file data/tcga_hdf5/methyTCGAfilt_assays.h5 \
--hdf5_se_file data/tcga_hdf5/methyTCGAfilt_se.rds --name methyTCGAfilt --step preprocess2D -profile docker -resume

## Random search CNN 2D
nextflow run workflows/train_model.nf --hdf5_file results/preprocess/methyTCGAfilt/2021-10-07/assay_reshaped.h5 \
--name CNN_autoencod --step random --network conf/network_params/cnn_network_autoencod_1.py \
--config conf/network_params/random_search_cnnautoencodv1_1.py -profile docker -resume

## Train CNN 2D
nextflow run workflows/train_model.nf --hdf5_file results/preprocess/methyTCGAfilt/2021-10-07/assay_reshaped.h5 \
--name CNN_autoencod --step train --network conf/network_params/cnn_network_autoencod_1.py \
--network_params conf/network_params/params_cnn_network_autoencod__v1_1.py -profile docker -resume

## Prepare TCGA data CNN 1D (test)
nextflow run workflows/train_model.nf --hdf5_file data/tcga_hdf5/testTCGA_assays.h5 \
--hdf5_se_file data/tcga_hdf5/testTCGA_se.rds --name testTCGA --step preprocessCNN -profile docker  -resume

## Run Random Search CNN 1D
nextflow run workflows/train_model.nf --hdf5_file results/preprocess/testTCGA/2021-10-26/assay_reshaped.h5 \
--name CNN_autoencod1D --step random --network conf/network_params/cnn1D_network_autoencod_1.py \
--config conf/network_params/random_search_cnn1D_autoencodv1_1.py -resume -profile docker --test_prop 0.8

nextflow run workflows/train_model.nf --hdf5_file results/preprocess/testTCGA/2021-10-26/assay_reshaped.h5 \
--name CNN_autoencod --params_name v1 --step train --network conf/network_params/cnn1D_network_autoencod_1.py \
--network_params conf/network_params/params_cnn1D_network_autoencod_v1_1.py  -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/preprocess/testTCGA/2021-10-26/assay_reshaped.h5 \
--name CNN_autoencod --params_name v1 --step transfer --layer "dense_1" --model results/CNN_autoencod/2021-10-28/model_trained/v1/CNN_autoencod/ \
--labels results/CNN_autoencod/2021-10-28/model_trained/v1/CNN_autoencod_labels.pb  -profile docker -resume

## Preprocess GSE97362
nextflow run workflows/train_model.nf --hdf5_file data/GSE97362/GSE97362_rawassays.h5 \
--hdf5_se_file data/GSE97362/GSE97362_rawse.rds --probes results/preprocess/methyTCGAfilt/2021-10-19/cpg_medians.Rdata \
  --name GSE97362 --step preprocessCNN_trans -profile docker -resume

nextflow run workflows/train_model.nf --hdf5_file results/preprocess/testTCGA/2021-10-26/assay_reshaped.h5 \
--name CNN_autoencod --params_name v1 --step autoencoder --model results/CNN_autoencod/2021-10-28/model_trained/v1/CNN_autoencod/ \
--labels results/CNN_autoencod/2021-10-28/model_trained/v1/CNN_autoencod_labels.pb  -profile docker -resume

## Test autoencoder
nextflow run workflows/train_model.nf --hdf5_file results/preprocess/testTCGA/2021-10-26/assay_reshaped.h5 --name dnn_autoencoder \
--params_name v1 --step autoencoder --network_params conf/network_params/params_dnn_autoencod_network1_v1.py  \
--network conf/network_params/dnn_autoencod_network1.py -resume -profile docker
