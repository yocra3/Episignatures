import pickle
import csv

from numpy import array
from numpy import argmax
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import numpy as np
from sklearn.model_selection import RandomizedSearchCV

with open('labels.txt','r') as file:
    project = file.read()
project = project.split('\n')[0:9813]
project = project[0:9813]

f = h5py.File('assays.h5', 'r')

## Prepare objects
meth_matrix = f['assay001']

## embedding labels
# integer encode
label_encoder = LabelEncoder()
integer_encoded = label_encoder.fit_transform(project)
# binary encode
onehot_encoder = OneHotEncoder(sparse=False)
integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
onehot_encoded = onehot_encoder.fit_transform(integer_encoded)

x_train = meth_matrix.reshape(meth_matrix.shape[0], meth_matrix.shape[1], 1)


## Model
def make_model(dense_layer_sizes, filters, kernel_size, stride):
    input_shape = (len(x_test[0]), 1)
    num_classes = len(y_train[0])
    
    model = Sequential()
    ## *********** First layer Conv
    model.add(Conv1D(filters, kernel_size = kernel_size, strides = min(stride, kernel_size), 
      input_shape = input_shape))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(2))
    model.output_shape

    ## ********* Classification layer
    model.add(Flatten())
    model.add(Dense(dense_layer_sizes, activation='relu'))
    model.add(Dense(num_classes, activation='softmax'))
    model.output_shape

    model.compile(loss = 'categorical_crossentropy',
                  optimizer = 'adam',
                  metrics = ['balanced_accuracy_score'])
    model.summary()
    return model
  

dense_size_candidates = [64, 128, 256, 512, 1024]
my_classifier = KerasClassifier(make_model, batch_size = 128)
validator = RandomizedSearchCV(my_classifier,
                         param_distributions = {'dense_layer_sizes': dense_size_candidates,
                                     # epochs is avail for tuning even when not
                                     # an argument to model building function
                                     'epochs': [25],
                                     'filters': [8, 16, 32, 64],
                                     'kernel_size': [(1, 1000)]},
                                     'stride': [(1, 1000)]
                         scoring = 'balanced_accuracy_score',
                         n_jobs = 1,
                         n_iter = 100, 
                         cv = 10,
                         random_state = 1)
search = validator.fit(x_train, onehot_encoded)

print('The parameters of the best model are: ')
print(search.best_params_)
# write it in a excel file
with open('results_runs50.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in validator.cv_results_.items():
       search.writerow([key, value])
# validator.best_estimator_ returns sklearn-wrapped version of best model.
# validator.best_estimator_.model returns the (unwrapped) keras model
best_model = search.best_estimator_.model
metric_names = best_model.metrics_names
metric_values = best_model.evaluate(x_test, y_test)
for metric, value in zip(metric_names, metric_values):
    print(metric, ': ', value)

pickle.dump( search, open( "model.p", "wb" ) )
