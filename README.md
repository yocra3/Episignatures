# Pathway activation inference using Machine Learning

This repository contains the code needed to perform the inference of pathway activation values using a Machine Learning framework. The repository contains the pipelines for two main operations:

1. Train a neural network to infere pathway activation
2. Compute the activation values for a new dataset

We provide a pre-trained network that returns the activation values for a selection of GO terms and KEGG pathways. 

## Requirements

- Nextflow
- Docker/Singularity
- R and Python (for data pre-processing)

## Pre-processing

Before using the pipeline for any of the steps, the data must be preprocessed. The input data is a h5 file containing a matrix with the gene expression values. 

In our pre-trained network, gene expression should be standardized. We provide a [script](scripts/create_HDF5_gexp.py) showing how to standardize the data and save it in the correct format.

When working with RNA-seq data, we recommend preprocessing the counts with the `DESeq2` Variant Stabilizing Transformation (vst) and then standardize the data. We show an example in this [script](scripts/download_gtex.R). For gene expression array data, expression data can be standardized without additional transformations. 

Finally, the input matrix should contain all the genes used to train the model. If some are not available, they should be filled with 0s, as shown in this [script](scripts/downloadGSE169038.R). 

## Training the model

We should use the nextflow pipeline to train the model. The following command shows how to call the pipeline:

```
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.6 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -profile docker
```

Pipeline arguments:

- `--hdf5_file`: h5 file with the gene expression data pre-processed. 
- `--name`: Name of the model. This name will be included in the folders and path names generated. We recoomend to used it to identify the data used to create the model. 
- `--params_name`: Name of the model configuration. This name will be included in the folders and path names generated. We recoomend to used it to identify the configuration of the model. 
- `--step`: For training the model, should be set to `train_biopathways`. 
- `--probes`: Names of the genes of `--hdf5_file`. They should be in the same order than in `--hdf5_file`. 
- `--autoencoder`: This value should always be set to `autoencoder` (still included for legacy issues).
- `--network_params`: Path to a python file including the parameters used to train the model (e.g. number of epochs for training, dropout percentage...).
- `--network`: Path to a python file including the neural network definition (see an [example](conf/network_params/dnn_gexp_autoencod_pathway_network3.py)).
- `--cpgmap`: Text file containing the match between genes and pathways.

Arguments starting by a single hyphen (-), belong to nextflow and can be found on its [documentation](https://www.nextflow.io/docs/latest/cli.html). 

This pipeline will output the results to 'results/name/params_name/model_trained' (If nextflow parameters are not changed). This folder will contain the following files:

- pathways_names.txt: names of the pathways included in the model.
- name/: folder with the trained model.
- name_history_model.pb: training information of the model. 
- name_training_evaluation.txt: training evaluation of the model (format txt).
- test_indices.csv: indices of the samples used for validation while training (0-based).

## Getting the activation values

Once we have the model trained, we can get the activation values with: 

```
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD -profile docker
```


Pipeline arguments:

- `--hdf5_file`: h5 file with the gene expression data pre-processed. 
- `--name`: Name of the dataset. This name will be included in the folders and path names generated. We recoomend to used it to identify the data used to create the model. 
- `--params_name`: Name of the model configuration. This name will be included in the folders and path names generated. We recoomend to used it to identify the configuration of the model. 
- `--step`: For training the model, should be set to `features`. 
- `--model`: Path to the model trained in the previous step.

This pipeline will output the results to 'results/name/params_name/model_features' (If nextflow parameters are not changed). This folder will only contain tsv files corresponding to the activation values of the intermediate layers. If our pre-trained network is used, only the file `prune_low_magnitude_dense.tsv` is generated. For other models with additional layers, other files can be generated. 

