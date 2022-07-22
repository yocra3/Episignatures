docker run -it -v /home/SHARED/PROJECTS/Episignatures:/home/SHARED/PROJECTS/Episignatures -w "$PWD" yocra3/episignatures_python:1.4  /bin/bash

python

#'#################################################################################
#'#################################################################################
#'  Retrain GTEx for TCGA
#'#################################################################################
#'#################################################################################

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

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model

## Train
A = open("results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/train.pb", 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

## Test
A = open("results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/test.pb", 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

## Model
model = tf.keras.models.load_model('results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding')
history = model.fit(x_train, y_train,
  batch_size = 128,
  epochs = 50,
  verbose = 1, validation_data = (x_test, y_test), workers = 10)
model.save( "results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/GTEX_transfer" )
