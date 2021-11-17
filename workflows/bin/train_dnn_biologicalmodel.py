#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Train a neural network defined in Keras
#'#################################################################################
#'#################################################################################

import pickle
import csv
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import h5py
import tensorflow as tf

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import concatenate, Dense, Dropout, Activation, Flatten, Input
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.regularizers import l2, l1

gpus = tf.config.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)


sys.path.append('./')
import params

name = sys.argv[1]
jobs = int(sys.argv[2])

f = h5py.File('assay_reshaped.h5', 'r')
proj_labels = f['label'].attrs['labels']
f.close()

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

with open('inputcpgs.txt','r') as file:
    cpgs = file.read()
cpgs = cpgs.split('\n')[0:-1]

## Select map rows with cpgs in input
with open('cpgs_map.txt','r') as file:
    cpg_map = pd.read_csv(file, delimiter = ' ', index_col = False)

cpg_map_filt = cpg_map[cpg_map['cpgs'].isin(cpgs)]
gene_count = cpg_map_filt['gene'].value_counts()
small_genes = gene_count.index[gene_count < 3]
mask = cpg_map_filt['gene'].isin(small_genes.to_list())

cpg_map_filt.loc[cpg_map_filt.index[mask], 'gene'] = "Intergenic"
gene_count2 = cpg_map_filt['gene'].value_counts()
cpg_map_filt['idx'] = [cpgs.index(cpg) for cpg in cpg_map_filt['cpgs']]

num_classes = len(y_train[0])

## Create list of inputs - gene
genes = gene_count2.index
genes = genes[genes != "Intergenic"]
gene_train = []
gene_test = []
gene_names = []
cpgs_input = []
for gene in genes:
  subset = cpg_map_filt[cpg_map_filt['gene'] == gene]
  selcpgs = subset['idx'].to_list()
  gene_train.append(x_train[:, selcpgs])
  gene_test.append(x_test[:, selcpgs])
  gene_names.append(gene)
  cpgs_input.append(subset['cpgs'].to_list())
  
## Create list of inputs - chr
if params.chr_genes:
  cpg_map_filt = cpg_map_filt[cpg_map_filt['gene'] == "Intergenic"]
chr_train = []
chr_test = []
chr_names = []
for chrom in cpg_map_filt['chromosome'].unique():
  subset = cpg_map_filt[cpg_map_filt['chromosome'] == chrom]
  selcpgs = subset['idx'].to_list()
  chr_train.append(x_train[:, selcpgs])
  chr_test.append(x_test[:, selcpgs])
  chr_names.append(chrom)
  cpgs_input.append(subset['cpgs'].to_list())

print("Creating gene model")
## Define model
gene_layer = []
for i in range(len(genes)):
  inputN = Input(shape=(gene_train[i].shape[1],))
  x = Dense(params.gene_neurons, activation = "relu", name = gene_names[i])(inputN)
  x = Model(inputs = inputN, outputs = x)
  gene_layer.append(x)

print("Creating chromosome model")
chr_layer = []
for i in range(len(chr_names)):
  inputN = Input(shape=(chr_train[i].shape[1],))
  x = Dense(params.chr_neurons, activation = "relu", name = chr_names[i])(inputN)
  x = Model(inputs = inputN, outputs = x)
  chr_layer.append(x)

layergene = [x.output for x in gene_layer]
layerchr = [x.output for x in chr_layer]
layer1 = concatenate(layergene + layerchr)

z = Dense(params.dense_layer_sizes1, activation="relu")(layer1)
z = Dense(params.dense_layer_sizes2, activation="relu")(z)
out = Dense(num_classes, activation='softmax')(z)

inp_gene = [x.input for x in gene_layer]
inp_chr = [x.input for x in chr_layer]
model = Model(inputs = inp_gene + inp_chr, outputs = out)

# Train model ####
opt = Adam(learning_rate = params.alpha)
model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
callbacks = [EarlyStopping(monitor = 'val_loss', patience = 5, verbose = 1)]

model.summary()
history = model.fit(gene_train + chr_train, y_train,
  batch_size = params.batch_size,
  epochs = params.epochs,
  verbose = 1, callbacks = callbacks, validation_data = (gene_test + chr_test, y_test), workers = jobs)

history_dict = history.history
pickle.dump( history_dict, open( name + "_history_model.pb", "wb" ), protocol = 4 )
pickle.dump( proj_labels, open( name + "_labels.pb", "wb" ), protocol = 4 )
model.save( name )

x_test = gene_test + chr_test
pickle.dump( [x_test, y_test], open( "test_list.pb", "wb" ), protocol = 4 )

pickle.dump( cpgs_input, open( "input_cpgs.pb", "wb" ), protocol = 4 )
