##################################################################################
##################################################################################
# DNN gexp network v3
##################################################################################
##################################################################################

from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import regularizers



def model_generator_train(x_train, y_train, params):
  num_classes = len(y_train[0])

  model = Sequential()
  model.add(Dense(params.dense_layer_sizes1, kernel_regularizer=regularizers.l1_l2(l1=params.l1, l2=params.l2), input_dim = x_train.shape[1], activation='relu'))
  model.add(Dense(params.dense_layer_sizes2, kernel_regularizer=regularizers.l1_l2(l1=params.l1, l2=params.l2), activation='relu'))
  model.add(Dense(num_classes, activation='softmax'))
  opt = Adam(learning_rate = params.alpha)
  model.compile(loss = 'categorical_crossentropy',
    optimizer = opt,
    metrics = ['accuracy'])
  return model
