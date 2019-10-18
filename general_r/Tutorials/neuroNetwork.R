#Demonstration of neuro networks using iris dataset

#ref: http://www.di.fc.ul.pt/~jpn/r/neuralnets/neuralnets.html
#The following code seems to work but I haven't fully understand the rationale
attach(iris)

library(caret)
library(dplyr)
library(e1071)
library(neuralnet)

set.seed(101)
in_train=createDataPartition(Species, p = 3/4)[[1]]
train_iris=iris[in_train,]
test_iris=iris[-in_train,]



modelFit <- train(Species ~ ., data=train_iris, method='nnet', trace = FALSE)
confusionMatrix(test_iris$Species,predict(modelFit,test_iris))


#manual calculation of accuracy
mean(as.character(predict(modelFit,test_iris))==as.character(test_iris$Species)) #0.972
