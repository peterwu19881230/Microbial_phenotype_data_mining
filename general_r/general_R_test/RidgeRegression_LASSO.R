#Ridge regression and the Lasso
#Ref: https://www.r-bloggers.com/ridge-regression-and-the-lasso/

?swiss
str(swiss) #47 obs. of  6 variables 

x=model.matrix(Fertility~.,swiss)[,-1]
y=swiss$Fertility

lambda=10^seq(10,-2,length=100)

library(glmnet)

set.seed(101)
train=sample(1:nrow(x),nrow(x)/2) #split to 50% train and 50% test
length(train)

test=-train
length(test)
ytest=y[test]

#OLS
swisslm=lm(Fertility~.,data=swiss)
coef(swisslm)

#ridge
ridge.mod=glmnet(x,y,alpha=0,lambda=lambda)
predict(ridge.mod,s=0,exact=T,type='coefficients')[1:6,] #this part that I copied don't run. No idea why.

swisslm = lm(Fertility~., data = swiss, subset = train)
ridge.mod = glmnet(x[train,], y[train], alpha = 0, lambda = lambda)
##find the best lambda from our list via cross-validation
cv.out = cv.glmnet(x[train,], y[train], alpha = 0)
## Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations
## per fold
bestlam = cv.out$lambda.min

#laso
lasso.mod=glmnet(x[train,],y[train],alpha=1,lambda=lambda)
lasso.pred=predict(lasso.mod,s=bestlam,newx=x[test,])

result=cbind(lasso.pred,ytest)
plot(result[,1],result[,2])
cor(result[,1],result[,2]) #0.76
#I still think I don't understand lasso completely


