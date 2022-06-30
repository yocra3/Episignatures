#'#################################################################################
#'#################################################################################
#' Check initial part of training with hsa00430
#'#################################################################################
#'#################################################################################

docker run -it -v /home/SHARED/PROJECTS/Episignatures:/home/SHARED/PROJECTS/Episignatures -w "$PWD" yocra3/episignatures_python:1.4  /bin/bash

python


import pickle
import csv
import os
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import h5py
import tensorflow as tf
import time

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input, Reshape
import tensorflow_model_optimization as tfmot

jobs = 10

A = open('results/GTEx_coding/paths_filt2_full_v3.6/model_trained/train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('results/GTEx_coding/paths_filt2_full_v3.6/model_trained/test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()


# Load genes
with open('results/GTEx_coding/input_genes.txt','r') as file:
    genes = file.read()
genes = genes.split('\n')[0:-1]

## Select map rows with cpgs in input
with open('results/GTEx_coding/go_kegg_final_gene_map.tsv','r') as file:
    gene_map = pd.read_csv(file, delimiter = '\t', index_col = False)


gene_map_filt = gene_map[gene_map['Symbol'].isin(genes)]
pathway_count = gene_map_filt['PathwayID'].value_counts()

gene_dict = dict(zip(genes, range(len(genes))))
gene_map_filt['idx'] = [gene_dict[gene] for gene in gene_map_filt['Symbol']]

## Create list of inputs - pathway
pathway_dict = {}
for g, idx in zip(gene_map_filt['PathwayID'].values, gene_map_filt['idx'].values):
  if g not in pathway_dict:
    pathway_dict[g] = [idx]
  else:
    pathway_dict[g].append(idx)

## Define gene mask ####
ngenes = len(genes)
npaths = len(pathway_dict.keys())

mask_d  = np.zeros((ngenes, npaths), float)

pcol = 0

pathways = pathway_dict.keys()
for path in pathways:
  mask_d[pathway_dict[path], pcol] = 1
  pcol += 1


path_train = x_train[:, pathway_dict['path:hsa00430']]
path_test = x_test[:, pathway_dict['path:hsa00430']]


opt = Adam(learning_rate = 0.01)
#
minimodel = Sequential()
minimodel.add(Dense(1, input_dim = path_train.shape[1], activation = 'elu'))
minimodel.add(Dense(path_train.shape[1]))
minimodel.compile(loss='mse',
    optimizer = opt,
    metrics = ['mse'])
minimodel.fit(path_train, path_train,
    batch_size = 128,
    epochs = 25,
    verbose = 1, validation_data = (path_test,
    path_test), workers = jobs)

ini_w = minimodel.get_weights()[0]

for i in range(9):
    minimodel = Sequential()
    minimodel.add(Dense(1, input_dim = path_train.shape[1], activation = 'elu'))
    minimodel.add(Dense(path_train.shape[1]))
    minimodel.compile(loss='mse',
        optimizer = opt,
        metrics = ['mse'])
    minimodel.fit(path_train, path_train,
        batch_size = 128,
        epochs = 25,
        verbose = 1, validation_data = (path_test,
        path_test), workers = jobs)
    ini_w = np.concatenate((ini_w, minimodel.get_weights()[0]), axis = 1)

f = h5py.File('results/GTEx_coding/paths_filt2_full_v3.6/hsa00430_ini_weights.h5', 'w')
dataset_input = f.create_dataset('weights_paths', (ini_w.shape[0], ini_w.shape[1]))
dataset_input[...] = ini_w
f.close()
