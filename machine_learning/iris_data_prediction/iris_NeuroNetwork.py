import numpy as np
from sklearn import datasets
from keras.layers import Dense #I did this in terminal first: 1. conda install keras 
from keras.models import Sequential
from keras.utils import to_categorical
import tensorflow as tf
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split


#For reproducible result I need the following 4 lines (Ref: https://machinelearningmastery.com/reproducible-results-neural-networks-keras/)
#================================
from numpy.random import seed
seed(101)
from tensorflow import set_random_seed
set_random_seed(102)
#================================



iris = datasets.load_iris()
X = iris.data
y = iris.target

#create traning set and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state=42,stratify=y) #not sure if stratify=y works when there are more than 2 classes for the labels


#instantiate the model + fit the data


##The following code is a combination of Ref1 and Ref2
##Ref 1 : https://campus.datacamp.com/courses/deep-learning-in-python/building-deep-learning-models-with-keras?ex=9
##Ref 2 : https://www.kaggle.com/akashsri99/deep-learning-iris-dataset-keras
model = Sequential()
model.add(Dense(10,activation='relu',input_shape=(X.shape[1],))) #X.shape[1] should equal the number of features
model.add(Dense(10,activation='relu'))
model.add(Dense(3,activation='softmax')) #the number of the output layer should equal the number of unique outcomes (response variables)

model.compile(optimizer='sgd',loss='categorical_crossentropy',metrics=['accuracy'])



model.fit(X_train, to_categorical(y_train),epochs=100)


y_pred = model.predict(X_test)
y_pred = np.argmax(y_pred,axis=1)


acc = accuracy_score(y_test, y_pred)

tf.Print(acc) 
