#! /opt/conda/envs/episignatures/bin/python

#'#################################################################################
#'#################################################################################
#' Export tcga model results
#'#################################################################################
#'#################################################################################


import pickle
import csv
import numpy as np
import sys
import pandas as pd
import h5py

from tensorflow.keras.models import Sequential
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import Sequential, Model
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint
from sklearn.metrics import confusion_matrix, classification_report


f = h5py.File('assay_reshaped.h5', 'r')
proj_labels = f['label'].attrs['labels']
f.close()


A = open('history_model.pb', 'rb')
history = pickle.load(A)
A.close()

A = open('model.pb', 'rb')
model = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

with open('training_evaluation.tsv', 'w') as csv_file:
    writer = csv.writer(csv_file, delimiter = '\t', escapechar=' ', quoting = csv.QUOTE_NONE)
    for key, value in history.items():
       writer.writerow([key + '\t' + '\t'.join([str(item) for item in value ])])


Y_pred = model.predict(x_test)
y_pred = np.argmax(Y_pred, axis=1)
y_class = np.argmax(y_test, axis=1)

cm = confusion_matrix(y_class, y_pred)
df = pd.DataFrame(cm, columns = proj_labels)
df.to_csv('confussionMatrix.tsv',  sep = "\t", index = False)

cr = classification_report(y_class, y_pred, target_names=proj_labels)
text_file = open("classificationReport.txt", "w")
text_file.write(cr)
text_file.close()
