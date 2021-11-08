#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#' Run random search on tcga data
#'#################################################################################
#'#################################################################################

import h5py
import pickle
import sys

from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

prop = float(sys.argv[1])


f = h5py.File('./assays.h5', 'r')
labels = f['label'][...]
methy = f['methy'][...]
f.close()

## embedding labels
# binary encode
onehot_encoder = OneHotEncoder(sparse=False)
integer_encoded = labels.reshape(len(labels), 1)
onehot_labels = onehot_encoder.fit_transform(integer_encoded)

x_train, x_test, y_train, y_test = train_test_split(methy, onehot_labels,
                                                    stratify = onehot_labels,
                                                    test_size = prop, random_state = 42)
pickle.dump( [x_train, y_train], open( "train.pb", "wb" ), protocol = 4 )
pickle.dump( [x_test, y_test], open( "test.pb", "wb" ), protocol = 4 )
