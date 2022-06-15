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

## filt2 - V6.2
model = tf.keras.models.load_model('results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD')
w = model.get_weights()

model2 = Sequential()
model2.add(Dense(w[0].shape[1], input_dim = w[0].shape[0], activation = 'elu'))
model2.add(Dense(w[5].shape[1]))

new_w = [w[i] for i in [0, 1, 5, 6]]
new_w[0] = w[0] * w[2]
new_w[2] = w[5]

model2.set_weights(new_w)

f = h5py.File('results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/model_weights.h5', 'w')
dataset_input = f.create_dataset('weights_paths', (new_w[0].shape[0], new_w[0].shape[1]))
dataset_input[...] = new_w[0]
f.close()

#
# def relprop_dense(self, x, w, r):
#     w_pos = tf.maximum(w, 0.0)
#     z = tf.matmul(x, w_pos) + self.epsilon
#     s = r / z
#     c = tf.matmul(s, tf.transpose(w_pos))
#     return c * x
#


## Adapted FROM github.com/KaiFabi/LayerwiseRelevancePropagation/
class RelevancePropagation(object):
    """Very basic implementation of the layer-wise relevance propagation algorithm.
    """
    def __init__(self, model, epsilon, rule):
        self.epsilon = epsilon
        self.rule = rule
        # Load model
        input_shape = model.input_shape
        # weights = model.get_weights()
        self.model = model
        # Extract model's weights
        self.weights = {weight.name.split('/')[0]: weight for weight in self.model.trainable_weights
                        if 'bias' not in weight.name}
        # # Extract activation layers
        self.activations = [layer.output for layer in self.model.layers]
        self.activations = self.activations[::-1]
        # self.activations.append(self.model.input)
        # Extract the model's layers name
        self.layer_names = [layer.name for layer in self.model.layers]
        # print(self.layer_names)
        self.layer_names = self.layer_names[::-1]
        self.layer_names = self.layer_names[:-1]
        print(self.layer_names)
    def run(self, input, index):
        # Build relevance graph
        self.relevance = self.relevance_propagation(index)
        f = K.function(self.model.input, self.relevance)
        # image = preprocess_input(input)
        # image = tf.expand_dims(input, axis=0)
        relevance_scores = f(input)
        # relevance_scores = self.postprocess(relevance_scores)
        return np.transpose(relevance_scores)
    def relevance_propagation(self, index):
        """Builds graph for relevance propagation."""
        relevance =  self.model.output[:, index:(index + 1)]
        for i, layer_name in enumerate(self.layer_names):
            if layer_name == self.layer_names[0]:
                w = self.weights[layer_name][:, index:(index + 1)]
            else:
                w = self.weights[layer_name]
            relevance = self.relprop_dense(self.activations[i+1], w, relevance)
        return relevance
    def relprop_dense(self, x, w, r):
        if self.rule == "z_plus":
            w_pos = tf.maximum(w, 0.0)
            z = tf.matmul(x, w_pos) + self.epsilon
            s = r / z
            c = tf.matmul(s, tf.transpose(w_pos))
            return c * x
        else:
            raise Exception("Error: rule for dense layer not implemented.")
    @staticmethod
    def rescale(x):
        """Rescales relevance scores of a batch of relevance maps between 0 and 1
        Args:
            x: RGB or grayscale relevance maps with dimensions (N, W, H, C) or (N, W, H), respectively.
        Returns:
            Rescaled relevance maps of same dimensions as input
        """
        return np.exp(x)/(np.sum(np.exp(x), axis = 1))
        # x_min = np.min(x, axis= 1)
        # x_max = np.max(x, axis= 1)
        # return (x - x_min).astype("float64") / (x_max - x_min).astype("float64")
    def postprocess(self, x):
        """Postprocesses batch of feature relevance scores (relevance_maps).
        Args:
            x: array with dimension (N, W, H, C)
        Returns:
            x: array with dimensions (N, W, H, C) or (N, W, H) depending on if grayscale or not
        """
        x = self.rescale(x)
        return x

f = h5py.File('results/HELIX/helix_array_reshaped_standardized.h5', 'r')
ori = f['methy'][...]
f.close()

## Define models by parts to reduce size of full matrix
model_pathways = Model(inputs=model2.input, outputs=model2.layers[0].output)

f = h5py.File('results/HELIX/comb_paths3_v3.6/relevance/relevance_genes_GO0032688.h5', 'w')
dataset_input = f.create_dataset('GO0032688', (ori.shape[0], ori.shape[1]))
lrp_pathways = RelevancePropagation(model_pathways, 1.0e-9, "z_plus")
dataset_input[...] = lrp_pathways.run(ori, 681).transpose()
f.close()


# f = h5py.File('results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/relevance/relevance_genes_pathways.h5', 'w')
# for i in range(model_pathways.output.shape[1]):
#     print(i)
#     dataset_input = f.create_dataset('patwhays_' + str(i), (ori.shape[0], ori.shape[1]))
#     model_pathways = Model(inputs=model2.input, outputs=model2.layers[1].output)
#     lrp_pathways = RelevancePropagation(model_pathways, 1.0e-9, "z_plus")
#     dataset_input[...] = lrp_pathways.run(ori, i).transpose()
#     del lrp_pathways
#     del dataset_input
#     del model_pathways
#     gc.collect()
#     K.clear_session()
# # f.close()


# f = h5py.File('results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/relevance/relevance_genes_pathways.h5', 'w')
# for i in range(model_pathways.output.shape[1]):
#     print(i)
#     dataset_input = f.create_dataset('patwhays_' + str(i), (ori.shape[0], ori.shape[1]))
#     model_pathways = Model(inputs=model2.input, outputs=model2.layers[1].output)
#     lrp_pathways = RelevancePropagation(model_pathways, 1.0e-9, "z_plus")
#     dataset_input[...] = lrp_pathways.run(ori, i).transpose()
#     del lrp_pathways
#     del dataset_input
#     del model_pathways
#     gc.collect()
#     K.clear_session()
# # f.close()
#
# ori_paths = model_pathways.predict(ori)
#
# f = h5py.File('results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/relevance/relevance_pathways_out.h5', 'w')
# for i in range(ori.shape[1]):
#     print(i)
#     dataset_input = f.create_dataset('out_' + str(i), (ori_paths.shape[0], ori_paths.shape[1]))
#     model_out = Model(inputs=model2.layers[1].output, outputs=model2.output)
#     lrp_out = RelevancePropagation(model_out, 1.0e-9, "z_plus")
#     dataset_input[...] =  lrp_out.run(ori_paths, i).transpose()
#     del lrp_out
#     del dataset_input
#     del model_out
#     gc.collect()
#     K.clear_session()
# f.close()
#
#
