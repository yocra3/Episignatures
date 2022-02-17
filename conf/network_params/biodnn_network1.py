##################################################################################
##################################################################################
# BIODNN network 1
##################################################################################
##################################################################################

import sys
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input, Reshape
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import regularizers
sys.path.append("/usr/GenNet/")
from GenNet_utils.LocallyDirectedConnected_tf2 import LocallyDirected1D


def model_generator_train(num_classes, gene_mask, params):

  inputsize = gene_mask.shape[0]
  inputN = Input(shape=(inputsize, ))
  Input_layer = Reshape(input_shape = (inputsize,), target_shape = (inputsize, 1))(inputN)
  Gene_layer = LocallyDirected1D(mask = gene_mask, filters = 1, input_shape=(inputsize, 1))(Input_layer)
  Gene_layer = Activation("relu")(Gene_layer)
  layer1 = Flatten()(Gene_layer)

  if params.dropout:
    layer1 = Dropout(params.dropout)(layer1)

  if params.regl1l2:
    z = Dense(params.dense_layer_sizes1, activation="relu", activity_regularizer=regularizers.l1_l2(l1=0.01, l2=0.01))(layer1)
  else:
    z = Dense(params.dense_layer_sizes1, activation="relu")(layer1)

  if params.dropout:
    z = Dropout(params.dropout)(z)

  if params.regl1l2:
    z = Dense(params.dense_layer_sizes2, activation="relu", activity_regularizer=regularizers.l1_l2(l1=0.01, l2=0.01))(z)
  else:
    z = Dense(params.dense_layer_sizes2, activation="relu")(z)

  out = Dense(num_classes, activation='softmax')(z)
  model = Model(inputs = inputN, outputs = out)
  return model



# def model_generator(x_train, y_train):
#     ## Define function to be passed
#     def model_function(dense_layer_sizes1, dense_layer_sizes2, dense_layer_sizes3, alpha):
#       num_classes = len(y_train[0])
#
#       model = Sequential()
#       ## *********** First layer Conv
#       model.add(Dense(dense_layer_sizes1, input_dim = x_train.shape[1], activation='relu'))
#       model.add(Dense(dense_layer_sizes2, activation='relu'))
#       model.add(Dense(dense_layer_sizes3, activation='relu'))
#       model.add(Dense(num_classes, activation='softmax'))
#       opt = Adam(learning_rate = alpha)
#       model.compile(loss = 'categorical_crossentropy',
#                     optimizer = opt,
#                     metrics = ['accuracy'])
#       model.summary()
#       return model
#     return model_function
