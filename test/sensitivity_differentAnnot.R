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
plot_metrics=function(samples,alpha,size,four_annot_list,metric,x_lab="high |PCC| -- ranked pairs -- low |PCC|",y){
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[1]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[1])),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[2]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[2])),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[3]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[3])),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[4]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[4])),size=size,alpha=alpha)+
    #geom_line(data=data.frame(no=samples,sensitivity=random_df$sensitivity[samples]),aes(no,sensitivity,color="sensitivity"),size=size,linetype = "dashed")+
    
    
    labs(x = x_lab,y=y,aesthetic='custom text',color="Metrics")+
    scale_x_continuous(expand = c(0,0),labels=comma)
}



#========================================================================================================================




#pwy
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcc")]
Pwy=cumsum(df$Pwy[order(df$pcc)])
Pwy_confusionMatrix_metrics=confusionMatrix_metrics(Pwy)
Pwy_metrics=Pwy_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]


#pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("pcomplex","pcc")]
pcomplex=cumsum(df$pcomplex[order(df$pcc)])
pcomplex_confusionMatrix_metrics=confusionMatrix_metrics(pcomplex)
pcomplex_metrics=pcomplex_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]


#pwy+pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","pcc")]
TF=( (df$Pwy+df$pcomplex)==2 )
sum(TF) #230
PwyANDpcomplex=cumsum( TF[order(df$pcc)])
PwyANDpcomplex_confusionMatrix_metrics=confusionMatrix_metrics(PwyANDpcomplex)
PwyANDpcomplex_metrics=PwyANDpcomplex_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]


#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
TF=( rowSums(df[,-6])==5 )
sum(TF) #106
All_annotSet=cumsum( TF[order(df$pcc)])
All_annotSet_confusionMatrix_metrics=confusionMatrix_metrics(All_annotSet)
All_annotSet_metrics=All_annotSet_confusionMatrix_metrics[,c("sensitivity","specificity","precision","accuracy")]




four_annot_list=list("EcoCyc Pathway"=Pwy_metrics,
                           "EcoCyc Protein complex"=pcomplex_metrics,
                           "EcoCyc Pathway and Protein complex"=PwyANDpcomplex_metrics,
                           "All 5 annotation sets"=All_annotSet_metrics)


##Set no. of samples to be plotted
##samples=1:200
##samples=1:4000
samples=seq(from=0,to=7914231,by=1989); samples[1]=1

p1=plot_metrics(samples=samples,alpha=0.5,size=1,four_annot_list,metric="sensitivity",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="Sensitivity")
p2=plot_metrics(samples=samples,alpha=0.5,size=1,four_annot_list,metric="specificity",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,four_annot_list,metric="precision",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,four_annot_list,metric="accuracy",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="accuracy")


require(gridExtra)
p_all=grid.arrange(p1, p2,p3,p4, ncol=2)


# (!) The following doesn't p_all to be properly generated
#pdf(file="metrics_diffAnnot_first4000",width=12*2,height=8*2) 
#p_all
#dev.off()

#pdf(file="metrics_diffAnnot",width=12,height=8)
#p_all
#dev.off()
