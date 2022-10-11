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
import keras_tuner as kt

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input, Reshape
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.utils import class_weight
from keras.regularizers import L1L2

# Load and reshape data
f = h5py.File('results/PRAD_metanalysis/prad_train_sel_assays.h5', 'r')
train_matrix = f['assay001']
train_matrix = train_matrix[...].reshape(train_matrix.shape[0], train_matrix.shape[1], 1)
f.close()

labels_df = pd.read_table("results/PRAD_metanalysis/individuals_labels.txt", delimiter=" ", names = ['Status', 'Project'])
labels_df['Comb'] = labels_df['Status'].astype(str) + labels_df['Project'].astype(str)

## embedding labels
# binary encode
label_encoder = LabelEncoder()
onehot_encoder = OneHotEncoder(sparse=False)

comb_int = label_encoder.fit_transform(labels_df['Comb'])
comb_encoded = comb_int.reshape(len(comb_int), 1)
onehot_comb = onehot_encoder.fit_transform(comb_encoded)

status_int = label_encoder.fit_transform(labels_df['Status'])

index = range(len(comb_int))


x_train, x_test, y_train, y_test, index_train, index_test = train_test_split(train_matrix, status_int, index,
                                                    stratify = onehot_comb,
                                                    test_size = 0.20, random_state = 42)
pickle.dump( [x_train, x_test, y_train, y_test, index_train, index_test], open( "results/PRAD_metanalysis/training_objects.pb", "wb" ), protocol = 4 )

A = open('results/PRAD_metanalysis/training_objects.pb', 'rb')
[x_train, x_test, y_train, y_test, index_train, index_test] = pickle.load(A)
A.close()


# Calculate the weights for each class so that we can balance the data
weights = class_weight.compute_class_weight('balanced', classes = np.unique(y_train),
                                            y = y_train)

weight = {0: weights[0], 1: weights[1]}
train_df = labels_df.loc[index_train]
n_project = train_df.groupby(['Project']).size()
samp_weight = 1/n_project[train_df.Project]

model = Sequential()

# model.add(Dense(200, activation='relu'))
# model.add(Dense(50, activation='relu'))
model.add(Dense(64, kernel_regularizer=L1L2(l1 = 0, l2 = 1e-6), bias_regularizer=L1L2(l1 = 1e-7, l2 = 1e-6), input_dim = x_train.shape[1], activation='relu'))
# model.add(Dense(200, kernel_regularizer=L1L2(l1 = 1e-6, l2 = 1e-6), bias_regularizer=L1L2(l1 = 1e-6, l2 = 1e-6), activation='relu'))
model.add(Dense(10, kernel_regularizer=L1L2(l1 = 1e-6, l2 = 0), bias_regularizer=L1L2(l1 = 0, l2 = 1e-7), activation='relu'))
# model.add(Dense(52, kernel_regularizer=L1L2(l1 = 1e-5, l2 = 1e-6), bias_regularizer=L1L2(l1 = 1e-5, l2 = 1e-4), activation='relu'))
# model.add(Dense(25, kernel_regularizer=L1L2(l1 = 1e-5, l2 = 1e-5), bias_regularizer=L1L2(l1 = 1e-5, l2 = 1e-5), activation='relu'))

# # model.add(Dropout(0.4))
# # model.add(Dropout(0.4))
# model.add(Dense(50, kernel_regularizer=L1L2(l1 = 1e-3, l2 = 1e-3), bias_regularizer=L1L2(l1 = 1e-3, l2 = 1e-3), activation='relu'))
model.add(Dense(1, activation='sigmoid'))


# Add the class weights to the training
opt = Adam(learning_rate = 1e-2)
model.compile(loss = 'binary_crossentropy',
    optimizer = opt,
    metrics = ['accuracy'])
model.summary()

history = model.fit(x_train, y_train,
  batch_size = 128,
  class_weight = weight,
  # sample_weight = samp_weight.to_numpy(),
  epochs = 30,
  verbose = 1, validation_data = (x_test, y_test), workers = 10)


Y_pred = model.predict(x_test)
y_pred = [0 if i < 0.5 else 1 for i in  Y_pred[:, 0]]

confusion_matrix(y_test, y_pred)

test_df = labels_df.loc[index_test]
test_df['pred'] = Y_pred
test_df.groupby(['Status', 'Project']).describe()

def model_builder(hp):
  model = Sequential()
  hp_units1 = hp.Int('units1', min_value=16, max_value=64, step=4)
  hp_units2 = hp.Int('units2', min_value=4, max_value=16, step=2)
  # hp_units3 = hp.Int('units3', min_value=16, max_value=64, step=4)
  hp_l1_k1 = hp.Choice('l1_k1', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  hp_l2_k1 = hp.Choice('l2_k1', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  hp_l1_b1 = hp.Choice('l1_b1', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  hp_l2_b1 = hp.Choice('l2_b1', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  hp_l1_k2 = hp.Choice('l1_k2', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  hp_l2_k2 = hp.Choice('l2_k2', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  hp_l1_b2 = hp.Choice('l1_b2', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  hp_l2_b2 = hp.Choice('l2_b2', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  # hp_l1_k3 = hp.Choice('l1_k3', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  # hp_l2_k3 = hp.Choice('l2_k3', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  # hp_l1_b3 = hp.Choice('l1_b3', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  # hp_l2_b3 = hp.Choice('l2_b3', values=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0])
  model.add(Dense(hp_units1, kernel_regularizer=L1L2(l1 = hp_l1_k1, l2 = hp_l2_k1), bias_regularizer=L1L2(l1 = hp_l1_b1, l2 = hp_l2_b1), input_dim = x_train.shape[1], activation='relu'))
  model.add(Dense(hp_units2, kernel_regularizer=L1L2(l1 = hp_l1_k2, l2 = hp_l2_k2), bias_regularizer=L1L2(l1 = hp_l1_b2, l2 = hp_l2_b2), activation='relu'))
  # model.add(Dense(hp_units3, kernel_regularizer=L1L2(l1 = hp_l1_k3, l2 = hp_l2_k3), bias_regularizer=L1L2(l1 = hp_l1_b3, l2 = hp_l2_b3), activation='relu'))
  model.add(Dense(1, activation='sigmoid'))
  # Tune the learning rate for the optimizer
  # Choose an optimal value from 0.01, 0.001, or 0.0001
  hp_learning_rate = hp.Choice('learning_rate', values=[1e-2, 1e-3, 1e-4, 1e-5])
  model.compile(optimizer= Adam(learning_rate=hp_learning_rate),
                loss='binary_crossentropy',
                metrics=['accuracy'])
  return model

tuner = kt.Hyperband(model_builder,
                     objective='val_loss',
                     max_epochs = 30,
                     hyperband_iterations = 5,
                     factor=3,
                     directory='results/PRAD_metanalysis/',
                     project_name='PRAD_metanalysis_sel')
tuner.search(x_train, y_train,   class_weight = weight, validation_split=0.2)

# classification_report(y_test, y_pred, target_names=['High', 'Low'])

# history_dict = history.history
# pickle.dump( history_dict, open( name + "_history_model.pb", "wb" ), protocol = 4 )
# pickle.dump( proj_labels, open( name + "_labels.pb", "wb" ), protocol = 4 )
