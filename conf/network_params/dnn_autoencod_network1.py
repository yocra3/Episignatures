##################################################################################
##################################################################################
# DNN autoencoder network 1
##################################################################################
##################################################################################

from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.optimizers import Adam

def model_generator_train(x_train, y_train, params):

  model = Sequential()
  model.add(Dense(params.dense_layer_sizes1, input_dim = x_train.shape[1], activation='relu'))
  model.add(Dense(params.dense_layer_sizes2, activation='relu'))
  model.add(Dense(params.dense_layer_sizes3, activation='relu'))
  model.add(Dense(params.dense_layer_sizes2, activation='relu'))
  model.add(Dense(params.dense_layer_sizes1, activation='relu'))
  model.add(Dense(x_train.shape[1], activation='sigmoid'))

  return model
