#capture sensitivity, accruacy and coverage on y-axis
#(In Yandell et. al 2012, it says: Three commonly used measures of gene-finder performance are sensitivity, specificity and accuracy)
#Another ref I personally prefer: https://www.dataschool.io/simple-guide-to-confusion-matrix-terminology/ (but Specificity is defined differently than the above paper)
#wikipedia: https://en.wikipedia.org/wiki/Confusion_matrix

#3 metrics defined in Yandell et. al 2012:
##sensitivity(SN) = TP/(TP+FN) 
##specificity(SP) = TP/(TP+FP)  #different from Wikipedia
##accuracy = (SN+SP/2) #different from Wikipedia. Also it looks weird. Doesn't make sense to me. I thought it should be: (TP + TN)/total 


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
  
  data.frame(
  TP=cumSums,
  FP=rankings-cumSums,
  TN=total-rankings-(max(cumSums)-cumSums),
  FN=max(cumSums)-cumSums,
  
  sensitivity=cumSums/max(cumSums),  # = TP/(TP+FN)
  specificity = (total-rankings-(max(cumSums)-cumSums)) / (total-max(cumSums)), # = TN / (TN+FP)
  precision= cumSums / rankings, # = TP/(TP+FP)
  accuracy =(cumSums+total-rankings-max(cumSums)+cumSums)/total # = (TP + TN)/total     ##Equation looks complicated. Have to double check
  )
  
}
#========================================================================================================================

#Function to get the graph. 
#========================================================================================================================
plot_metrics=function(samples,alpha,size,metric_df,random_df,x_lab="high |PCC| -- ranked pairs -- low |PCC|",y=""){
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
    geom_line(data=data.frame(no=samples,sensitivity=metric_df$sensitivity[samples]),aes(no,sensitivity,color="sensitivity"),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,specificity=metric_df$specificity[samples]),aes(no,specificity,color="specificity"),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,precision=metric_df$precision[samples]),aes(no,precision,color="precision"),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,accuracy=metric_df$accuracy[samples]),aes(no,accuracy,color="accuracy"),size=size,alpha=alpha)+
    
    geom_line(data=data.frame(no=samples,sensitivity=random_df$sensitivity[samples]),aes(no,sensitivity,color="sensitivity"),size=size,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,specificity=random_df$specificity[samples]),aes(no,specificity,color="specificity"),size=size,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,precision=random_df$precision[samples]),aes(no,precision,color="precision"),size=size,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,accuracy=random_df$accuracy[samples]),aes(no,accuracy,color="accuracy"),size=size,linetype = "dashed")+
    
    labs(x = x_lab,y=y,aesthetic='custom text',color="Metrics")+
    scale_x_continuous(expand = c(0,0),labels=comma)
}
#========================================================================================================================

##Set no. of samples to be plotted
##samples=1:4000
samples=seq(from=0,to=7914231,by=1989); samples[1]=1


#pwy
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcc")]
Pwy=cumsum(df$Pwy[order(df$pcc)])
Pwy_confusionMatrix_metrics=confusionMatrix_metrics(Pwy)
Pwy_metrics=Pwy_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]

set.seed(102)
Pwy_randomOrder=cumsum(sample(df$Pwy))
Pwy_randomOrder_confusionMatrix_metrics=confusionMatrix_metrics(Pwy_randomOrder)
Pwy_randomOrder_metrics=Pwy_randomOrder_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]
  


plot_metrics(samples=samples,alpha=0.5,size=1,metric_df=Pwy_metrics,random_df=Pwy_randomOrder_metrics) #specificity overlaps with accuracy 


  
#pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("pcomplex","pcc")]
pcomplex=cumsum(df$pcomplex[order(df$pcc)])
pcomplex_confusionMatrix_metrics=confusionMatrix_metrics(pcomplex)
pcomplex_metrics=pcomplex_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]

set.seed(102)
pcomplex_randomOrder=cumsum(sample(df$pcomplex))
pcomplex_randomOrder_confusionMatrix_metrics=confusionMatrix_metrics(pcomplex_randomOrder)
pcomplex_randomOrder_metrics=pcomplex_randomOrder_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]


plot_metrics(samples=samples,alpha=0.5,size=1,metric_df=pcomplex_metrics,random_df=pcomplex_randomOrder_metrics) #specificity overlaps with accuracy 


#pwy+pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","pcc")]
TF=( (df$Pwy+df$pcomplex)==2 )
sum(TF) #230
PwyANDpcomplex=cumsum( TF[order(df$pcc)])
PwyANDpcomplex_confusionMatrix_metrics=confusionMatrix_metrics(PwyANDpcomplex)
PwyANDpcomplex_metrics=PwyANDpcomplex_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]


set.seed(102)
PwyPcomplex_randomOrder=cumsum(sample(TF))
PwyPcomplex_randomOrder_confusionMatrix_metrics=confusionMatrix_metrics(PwyPcomplex_randomOrder)
PwyPcomplex_randomOrder_metrics=PwyPcomplex_randomOrder_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]


plot_metrics(samples=samples,alpha=0.5,size=1,metric_df=PwyANDpcomplex_metrics,random_df=PwyPcomplex_randomOrder_metrics) #specificity overlaps with accuracy 


#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
TF=( rowSums(df[,-6])==5 )
sum(TF) #106
All_annotSet=cumsum( TF[order(df$pcc)])
All_annotSet_confusionMatrix_metrics=confusionMatrix_metrics(All_annotSet)
All_annotSet_metrics=All_annotSet_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]


set.seed(102)
All_annotSet_randomOrder=cumsum(sample(TF))
All_annotSet_randomOrder_confusionMatrix_metrics=confusionMatrix_metrics(All_annotSet_randomOrder)
All_annotSet_randomOrder_metrics=All_annotSet_randomOrder_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]


plot_metrics(samples=samples,alpha=0.5,size=1,metric_df=All_annotSet_metrics,random_df=All_annotSet_randomOrder_metrics) #specificity overlaps with accuracy 









