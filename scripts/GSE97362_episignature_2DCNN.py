#'#################################################################################
#'#################################################################################
#'  Test 2D CNN in gse97362 using subsets of features
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

from numpy import array, argmax
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.layers import Conv2D, MaxPooling2D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import confusion_matrix, classification_report
from keras.regularizers import l2

label_encoder = LabelEncoder()
onehot_encoder = OneHotEncoder(sparse=False)

## Load methylation
f = h5py.File("results/gse97362_episignature_validation/mini_conv_2D.h5", 'r')
cot = f['discov'][...]
x_train = cot.reshape(cot.shape[0], cot.shape[1], cot.shape[2], 1)
cot = f['val'][...]
x_test = cot.reshape(cot.shape[0], cot.shape[1], cot.shape[2], 1)
f.close()

# Load labels and ccnvert them to integers 
with open('results/gse97362_episignature_validation/discov_labels.txt','r') as file:
    label = file.read()
label = label.split('\n')[0:-1]
int1 = label_encoder.fit_transform(label)
int_dim = int1.reshape(len(int1), 1)
y_train = onehot_encoder.fit_transform(int_dim)

with open('results/gse97362_episignature_validation/val_labels.txt','r') as file:
    label = file.read()
label = label.split('\n')[0:-1]
int1 = label_encoder.fit_transform(label)
int_dim = int1.reshape(len(int1), 1)
y_test = onehot_encoder.fit_transform(int_dim)

input_shape = (x_train.shape[1], x_train.shape[2], 1)
num_classes = len(y_train[0])

model = Sequential()
model.add(Conv2D(4, kernel_size = (7, 4), 
    strides = 1,
    input_shape = input_shape))
model.add(Activation('relu'))
model.add(MaxPooling2D((1, 4)))
## ********* Classification layer
model.add(Flatten())
model.add(Dense(2048, activation='relu'))
model.add(Dense(128, activation='relu'))
model.add(Dense(num_classes, activation='softmax'))

opt = Adam(learning_rate = 1e-3)
model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
history = model.fit(x_train, y_train,
  batch_size = 32,
  epochs = 50,
  verbose = 1, validation_data = (x_test, y_test))

Y_pred = model.predict(x_test)
y_pred = np.argmax(Y_pred, axis=1)
y_class = np.argmax(y_test, axis=1)
confusion_matrix(y_class, y_pred)

f = h5py.File('./results/gse97362_episignature_validation/conv_2D_outputs.h5', 'w')
for i in range(3, len(model.layers) - 1):
  new_model = Model(inputs=model.input, outputs=model.layers[i].output)
  Y_pred = new_model.predict(x_test)
  res = f.create_dataset(model.layers[i].name, (Y_pred.shape[0], Y_pred.shape[1]))
  res[...] = Y_pred
f.close()


f = h5py.File('./results/gse97362_episignature_validation/conv_2D_outputs_train.h5', 'w')
for i in range(3, len(model.layers) - 1):
  new_model = Model(inputs=model.input, outputs=model.layers[i].output)
  Y_pred = new_model.predict(x_train)
  res = f.create_dataset(model.layers[i].name, (Y_pred.shape[0], Y_pred.shape[1]))
  res[...] = Y_pred
f.close()
