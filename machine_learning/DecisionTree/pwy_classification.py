#This script preps the required .RData: pwy_neuroNetwork_prep.R

#Load the .RData (based on: https://stackoverflow.com/questions/21288133/loading-rdata-files-into-python)


#Do this in terminal: pip install pyreadr

import numpy as np
import pandas as pd
import pyreadr
import os
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix

os.chdir("/Users/peterwu/Google Drive File Stream/Shared Drives/HuLab/Nichols_Data_mining/")
objects = pyreadr.read_r('Data/phenoForPythonMachineLearning.RData')
pwyNewPheno=objects["pwyNewPheno"]
pwyLabel=objects["pwyLabel"]

type(pwyNewPheno) #pandas.core.frame.DataFrame
type(pwyLabel) #pandas.core.frame.DataFrame

pwyNewPheno.shape #(1954, 324)
pwyLabel.shape #(2344, 1)


#create traning set and test set
X=pwyNewPheno
y=pwyLabel

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state=42,stratify=y) #not sure if stratify=y works when there are more than 2 classes for the labels


#instantiate the model + fit the data


##I want to try dicision tree first
dt = DecisionTreeClassifier(random_state=43)

dt.fit(X_train, y_train)



#predict the test data and calculate accuracy, sensitivity...etc

y_pred=dt.predict(X_test)


acc = accuracy_score(y_test, y_pred) #0.018 (completely failed)



#Use label matrix
#===============================================================================

#Do this in terminal: pip install pyreadr
import pyreadr
import os
exec(open("/Users/peterwu/Google Drive File Stream/Shared Drives/HuLab/Nichols_Data_mining/test/machine_learning/DecisionTreePipeline.py").read()) #import the function I defined

os.chdir("/Users/peterwu/Google Drive File Stream/Shared Drives/HuLab/Nichols_Data_mining/")


objects = pyreadr.read_r('Data/sourced/All_Data_NAimputed.RData')
NewPheno=objects["All_Data_NAimputed"]


objects = pyreadr.read_r('Data/pwyLabel_df.RData')
Label=objects["pwyLabel_df"]




type(NewPheno) #pandas.core.frame.DataFrame
type(Label) #pandas.core.frame.DataFrame

NewPheno.shape 
Label.shape 





metrics_df=DecisionTreePipeline(NewPheno,Label)
metrics_df.to_csv("test/machine_learning/pathway.csv")