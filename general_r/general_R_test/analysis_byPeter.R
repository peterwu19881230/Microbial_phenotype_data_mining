#Build a model to distinguish CDA and PDA by microarray data. Method: Laso (https://www.youtube.com/watch?v=fAPCaue8UKQ)


library(dplyr); library(caret); library(e1071); library(neuralnet); library(readxl); library(glmnet)

#Load and clean the data
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)

dat=read_excel(paste(currentDir,"Norm_NR.xlsx",sep="/"),
               col_types=rep(c("text","numeric"),times=c(3,35)),
               na="NA"
               ) #got warnings( because of NA values). Also, columns with NA are automatically turned into character vector


dat=dat[,-(1:3)]
str(dat)

dat=as.matrix(dat)
str(dat)

#impute NA values with mean of each column
Index=(1:35)[apply(dat,2,anyNA)] #indices for columns containing NA

for(i in Index){
  imputation=mean(as.numeric(dat[,i]),na.rm=T) #calculate column mean for imputation
  dat[,i][is.na(dat[,i])]=imputation #impute the mean
}

#check to see if all NAs are imputed:
apply(dat,2,anyNA) %>% sum #0 NA

#transpose so all the genes become variable, CDA, PDA become the outcomes
dat=t(dat)
str(dat)


#Model building

response=str_extract(rownames(dat),pattern="[A-Z]{3}")
predictor=ifelse(response=="PDA",1,0)



start.time=Sys.time()
cv.lasso = cv.glmnet(dat, predictor, family='binomial', alpha=1, parallel=TRUE, standardize=TRUE, type.measure='class')
end.time=Sys.time()
end.time-start.time #Time difference of 3.486802 secs


fit=glmnet(dat, predictor, family='binomial', alpha=1,lambda=cv.lasso$lambda.1se)


predict(fit,dat) #Weird result. I don't know what's going on. I also have double that maybe the no. of patients are not enough




