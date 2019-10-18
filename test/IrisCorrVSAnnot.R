#Try to use iris dataset + pcc to do Corr VS. Annot (Sensitivity + Specificity + Precision + Accuracy)



#x: ranked pairs , y: Sensitivity --- different combinatio of annotation sets on the same graph
#(From confusionMatrix_analysis.R I noticed that sensitivuty is probably the only thing that is worth plotting)



#4 metrics defined in Wikipedia:
##sensitivity = TP / (TP +FN)
##specificity = TN / (TN+FP)
##precision = TP/(TP+FP)
##accuracy = (TP + TN)/total


##suppose the cutoff pcc=rho

##TP=no. of pairs that: (a) have pcc >= rho (b) are co-annotated
##TN=no. of pairs that: (a) have pcc < rho (b) are not co-annotated
##FP=no. of pairs that: (a) have pcc >= rho (b) are not co-annotated
##FN=no. of pairs that: (a) have pcc < rho (b) are co-annotated



#function to obtain confusion matrix and some metrics
#========================================================================================================================
confusionMatrix_metrics=function(cumSums,rankings=seq_along(cumSums),total=7914231){ #cumSums would be a numeric vector
  
  df=data.frame(
    TP=cumSums,
    FP=rankings-cumSums,
    TN=total-rankings-(max(cumSums)-cumSums),
    FN=max(cumSums)-cumSums,
    
    sensitivity=cumSums/max(cumSums),  # = TP/(TP+FN)
    specificity = (total-rankings-(max(cumSums)-cumSums)) / (total-max(cumSums)), # = TN / (TN+FP)
    precision= cumSums / rankings, # = TP/(TP+FP)
    accuracy =(cumSums+total-rankings-max(cumSums)+cumSums)/total # = (TP + TN)/total     ##Equation looks complicated. Have to double check
    
  )
  
  randomCoAnnotation=rep(0,total); set.seed(102); randomCoAnnotation[sample(total,max(cumSums))]=1
  randomCumSums=cumsum(randomCoAnnotation)
  
  random_df=data.frame(
    random_TP=randomCumSums,
    random_FP=rankings-randomCumSums,
    random_TN=total-rankings-(max(randomCumSums)-randomCumSums),
    random_FN=max(randomCumSums)-randomCumSums,
    
    random_sensitivity=randomCumSums/max(randomCumSums),  # = TP/(TP+FN)
    random_specificity = (total-rankings-(max(randomCumSums)-randomCumSums)) / (total-max(randomCumSums)), # = TN / (TN+FP)
    random_precision= randomCumSums / rankings, # = TP/(TP+FP)
    random_accuracy =(randomCumSums+total-rankings-max(randomCumSums)+randomCumSums)/total # = (TP + TN)/total     ##Equation looks complicated. Have to double check
    
  )
  
  return(cbind(df,random_df))
  
}
#========================================================================================================================


#Function to get the graph. 
#========================================================================================================================
plot_metrics=function(samples,alpha,size,annot_df,metric,x_lab="high |PCC| -- ranked pairs -- low |PCC|",y,legend_title="Annotation set(s)"){
  
  random_metric=paste("random",metric,sep="_")
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
    geom_line(data=data.frame(no=samples,Metric=annot_df[[metric]][samples]),aes(no,Metric,color="Iris"),size=size,alpha=alpha)+
      
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=annot_df[[random_metric]][samples]),aes(no,Metric,color="Iris"),size=size,alpha=alpha+0.3,linetype = "dashed")+
     
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)
}

#========================================================================================================================


head(iris)
my_iris=iris

#center + scale the data
my_iris$Sepal.Length=scale(my_iris$Sepal.Length)
my_iris$Sepal.Width=scale(my_iris$Sepal.Width)
my_iris$Petal.Length=scale(my_iris$Petal.Length)
my_iris$Petal.Width=scale(my_iris$Petal.Width)
my_iris$Species=as.character(my_iris$Species)

dat=my_iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]
attribute_list=attr_list(rownames(my_iris),my_iris$Species)
result=dist_TF_cumsum_matirxOperation(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
cumSums=result$cumsum


iris_confusionMatrix_metrics=confusionMatrix_metrics(cumSums=cumSums,total=length(cumSums))

##samples=1:length(cumSums)
samples=1:1000
p1=plot_metrics(samples=samples,alpha=0.5,size=1,annot_df=iris_confusionMatrix_metrics,metric="sensitivity",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="Sensitivity")
p2=plot_metrics(samples=samples,alpha=0.5,size=1,annot_df=iris_confusionMatrix_metrics,metric="specificity",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,annot_df=iris_confusionMatrix_metrics,metric="precision",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,annot_df=iris_confusionMatrix_metrics,metric="accuracy",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="accuracy")


require(gridExtra)
p_all=grid.arrange(p1, p2,p3,p4, ncol=2) #When viewed in small window the dashed lines don't show properly


#from the low sensitivity compared to random control, pcc might not be the good way to group different iris species
