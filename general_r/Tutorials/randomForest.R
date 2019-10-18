#Demonstrate Random Forest Classification using Iris Data


#The following code seems to work but I haven't fully understand the rationale
#I coded the following based on the dreaded, unclear chapter of Coursera Data Science Lecture 8
attach(iris)

library(caret)
library(dplyr)
library(e1071)

set.seed(101)
in_train=createDataPartition(Species, p = 3/4)[[1]]
train_iris=iris[in_train,]
test_iris=iris[-in_train,]

modFit=train(Species~.,method="rf",prox=T,data=train_iris)
print(modFit)

confusionMatrix(Species[-in_train],predict(modFit,test_iris))


#manual calculation of accuracy
mean(as.character(predict(modFit,test_iris))==as.character(test_iris$Species)) #0.972 (same as using neuronetwork)
