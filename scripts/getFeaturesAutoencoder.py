docker run -v $PWD:$PWD -it  yocra3/episignatures_python:1.3

## marvin
srun --mem 256Gb --pty bash   
singularity exec -B /gpfs42/robbyfs/scratch/lab_laperez/cruizg/Episignatures -B "$PWD" /gpfs42/robbyfs/scratch/lab_laperez/cruizg/Episignatures/cache_singularity/yocra3-episignatures_python-1.3.img /bin/bash    

cd /gpfs42/robbyfs/scratch/lab_laperez/cruizg/Episignatures
python 

#'#################################################################################
#'#################################################################################
#'  Get features from autoencoder network for GSE97632
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

f = h5py.File('./results/preprocess/GSE97362/2021-10-29/assay_reshaped.h5', 'r')
methy = f['methy'][...]
f.close()

auto = tf.keras.models.load_model('results/CNN_autoencod/2021-11-04/autoencoder_trained/v1/CNN_autoencod_autoencoder') 
cnn = tf.keras.models.load_model('results/CNN_autoencod/2021-10-28/model_trained/v1/CNN_autoencod')

cnn_part = Model(inputs=cnn.input, outputs=cnn.layers[3].output)
methy_cnn = cnn_part.predict(methy)
methy_auto = auto.predict(methy_cnn)

model_d1 = Model(inputs=auto.input, outputs=auto.layers[0].output)
methy_dense1 = model_d1.predict(methy_cnn)
pd.DataFrame(methy_dense1).to_csv('results/GSE97632/2021-11-05/model_features/autoencoder_dense.tsv',  sep = "\t", index = False)

model_d2 = Model(inputs=auto.input, outputs=auto.layers[1].output)
methy_dense2 = model_d2.predict(methy_cnn)
pd.DataFrame(methy_dense2).to_csv('results/GSE97632/2021-11-05/model_features/autoencoder_dense1.tsv',  sep = "\t", index = False)

df = pd.DataFrame(methy_auto)
df.to_csv('results/GSE97632/2021-11-05/model_features/autoencoder_output.tsv',  sep = "\t", index = False)
