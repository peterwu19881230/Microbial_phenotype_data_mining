#Try to determine which subset of Nichols' can be used as gold standard to do statistical learning



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
  
  df=dfa.frame(
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
  
  random_df=dfa.frame(
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
    geom_line(data=data.frame(no=samples,Metric=annot_df[[metric]][samples]),aes(no,Metric,color="All 5 annotation sets"),size=size,alpha=alpha)+
      
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=annot_df[[random_metric]][samples]),aes(no,Metric,color="All 5 annotation sets"),size=size,alpha=alpha+0.3,linetype = "dashed")+
     
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)
}

#========================================================================================================================


names(id_allAttributes)

Subset=id_allAttributes[,c("ids","Pwy","pcomplex","regulator","operon","kegg_modules")] %>% unique
Subset=Subset[!apply(Subset,1,anyNA),]
dim(Subset)

idIn5AnnotationSets=Subset$ids %>% unique
length(idIn5AnnotationSets)


df=strain1strain2_allAnnotations_allDistances
df=df[ (as.character(df$strain1) %in% idIn5AnnotationSets )|
           (as.character(df$strain2) %in% idIn5AnnotationSets ) ,
         ]

dim(df)


#pwy + pwcomplex + KEGG modules + regulon + operon
TF=( rowSums(df[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5 )
sum(TF) #106

All_annotSet=cumsum( TF[order(df$pcc)])
All_annotSet_confusionMatrix_metrics=confusionMatrix_metrics(cumSums=All_annotSet,total=dim(df)[1])


##Set no. of samples to be plotted
##samples=1:200
##samples=1:4000 #I should use less than 4000 because dim(df)[1]<<7914231, but showing more doesn't hurt
##samples=1:993171 #too slow

p1=plot_metrics(samples=samples,alpha=0.5,size=1,annot_df=All_annotSet_confusionMatrix_metrics,metric="sensitivity",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="Sensitivity")
p2=plot_metrics(samples=samples,alpha=0.5,size=1,annot_df=All_annotSet_confusionMatrix_metrics,metric="specificity",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,annot_df=All_annotSet_confusionMatrix_metrics,metric="precision",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,annot_df=All_annotSet_confusionMatrix_metrics,metric="accuracy",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="accuracy")


require(gridExtra)
p_all=grid.arrange(p1, p2,p3,p4, ncol=2) #When viewed in small window the dashed lines don't show properly


#dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
#ggsave(file=paste(dir_of_workingScript,"subsetAsGoldStandard.pdf",sep="/"),plot=p_all,device="pdf",width=3*7,height=2*7)




