#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Extract features from bio DNN network
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
from tensorflow.keras.models import Model

A = open('input_list.pb', 'rb')
input_list =  pickle.load(A)
A.close()

model = tf.keras.models.load_model('model/' + os.listdir('model/')[0]) 

for i in range(len(model.layers) - 4, len(model.layers) - 1):
  new_model = Model(inputs=model.input, outputs=model.layers[i].output)
  Y_pred = new_model.predict(input_list)
  df = pd.DataFrame(Y_pred)
  df.to_csv(model.layers[i].name + '.tsv',  sep = "\t", index = False)


