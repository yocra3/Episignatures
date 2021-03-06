#'#################################################################################
#'#################################################################################
#' Convert a network from TF 2.7 to TF 2.2 for innvestigate
#'#################################################################################
#'#################################################################################

#' Save weights from original model
#'#################################################################################
docker run -it -v /home/SHARED/PROJECTS/Episignatures:/home/SHARED/PROJECTS/Episignatures -w "$PWD" yocra3/episignatures_python:1.4  /bin/bash

python


import pickle
import h5py
import csv
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import tensorflow as tf
import os
import tensorflow.keras.backend as K
from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense
import copy
import gc

## Main model GTEx
model = tf.keras.models.load_model('results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding')
w = model.get_weights()

new_w = [w[i] for i in [0, 1, 5, 6]]
new_w[0] = w[0] * w[2]
new_w[2] = w[5]

f = h5py.File('results/GTEx_coding/paths_filt2_full_v3.11/model_trained/model_weights.h5', 'w')
dataset_input = f.create_dataset('weights_paths', (new_w[0].shape[0], new_w[0].shape[1]))
dataset_input[...] = new_w[0]
f.close()

sub_mods = [chr(x) for x in range(97, 102)]

for sub in sub_mods:
    mod = tf.keras.models.load_model('results/GTEx_coding/paths_filt2_full_v3.11' + sub + '/model_trained/GTEx_coding')
    w_m = mod.get_weights()
    new_wm = w_m[0] * w_m[2]
    f = h5py.File('results/GTEx_coding/paths_filt2_full_v3.11' + sub + '/model_trained/model_weights.h5', 'w')
    dataset_input = f.create_dataset('weights_paths', (new_wm.shape[0], new_wm.shape[1]))
    dataset_input[...] = new_wm
    f.close()

model_pre = tf.keras.models.load_model('results/GTEx_coding/paths_filt2_pre_v3.8/model_trained/GTEx_coding')
w_pre = model_pre.get_weights()

new_wm_pre = w_pre[0] * w_pre[2]
f = h5py.File('results/GTEx_coding/paths_filt2_pre_v3.8/model_trained//model_weights.h5', 'w')
dataset_input = f.create_dataset('weights_paths', (new_wm_pre.shape[0], new_wm_pre.shape[1]))
dataset_input[...] = new_wm_pre
f.close()


for sub in sub_mods:
    mod = tf.keras.models.load_model('results/GTEx_coding/paths_filt2_pre_v3.8' + sub + '/model_trained/GTEx_coding')
    w_m = mod.get_weights()
    new_wm = w_m[0] * w_m[2]
    f = h5py.File('results/GTEx_coding/paths_filt2_pre_v3.8' + sub + '/model_trained/model_weights.h5', 'w')
    dataset_input = f.create_dataset('weights_paths', (new_wm.shape[0], new_wm.shape[1]))
    dataset_input[...] = new_wm
    f.close()

model_tcga = tf.keras.models.load_model('results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/TCGA_coding_all')
w_tcga = model_tcga.get_weights()

new_wm_tcga = w_tcga[0] * w_tcga[2]
f = h5py.File('results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/model_weights.h5', 'w')
dataset_input = f.create_dataset('weights_paths', (new_wm_tcga.shape[0], new_wm_tcga.shape[1]))
dataset_input[...] = new_wm_tcga
f.close()


model_tcga_pre = tf.keras.models.load_model('results/TCGA_coding_all/paths_filt2_pre_v3.8/model_trained/TCGA_coding_all')
w_tcga_pre = model_tcga_pre.get_weights()

new_wm_tcga_pre = w_tcga_pre[0] * w_tcga_pre[2]
f = h5py.File('results/TCGA_coding_all/paths_filt2_pre_v3.8/model_trained/model_weights.h5', 'w')
dataset_input = f.create_dataset('weights_paths', (new_wm_tcga_pre.shape[0], new_wm_tcga_pre.shape[1]))
dataset_input[...] = new_wm_tcga_pre
f.close()

model_tcga_retrain = tf.keras.models.load_model('results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/GTEX_transfer')
w_tcga_retrain = model_tcga_retrain.get_weights()

new_wm_tcga_retrain = w_tcga_retrain[0] * w_tcga_retrain[2]
f = h5py.File('results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/model_weights_retrain.h5', 'w')
dataset_input = f.create_dataset('weights_paths', (new_wm_tcga_retrain.shape[0], new_wm_tcga_retrain.shape[1]))
dataset_input[...] = new_wm_tcga_retrain
f.close()
