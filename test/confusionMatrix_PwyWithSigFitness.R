#Use conditions enriched by pathways to perform confusion matrix analysis



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
plot_metrics=function(samples,alpha,size,dist_list,metric,x_lab="Low distance -- ranked pairs -- High distance",y,legend_title="Types of data - distance"){
  
  random_metric=paste("random",metric,sep="_")
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
    geom_line(data=data.frame(no=samples,Metric=dist_list[[1]][[metric]][samples]),aes(no,Metric,color=names(dist_list[1])),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[2]][[metric]][samples]),aes(no,Metric,color=names(dist_list[2])),size=size,alpha=alpha)+
    
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=dist_list[[1]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[1])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[2]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[2])),size=size,alpha=alpha+0.3,linetype = "dashed")+
     
    
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)
}
#========================================================================================================================



df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcc")]
sum(df$Pwy) #7787
Original_pcc=cumsum( df$Pwy[order(df$pcc)])
Original_pcc_confusionMatrix_metrics=confusionMatrix_metrics(Original_pcc)



#Calculate confusion matrix and other metrices using the Pwys that have at least 1 sig phenotype
#================================================================================================
 



## Remove those pwys by assigning NA 
#------------------------------------------------------------------------------------
id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique

PwyLeft=id_allAttributes$Pwy[id_allAttributes$AnySigFitness & !is.na(id_allAttributes$Pwy)] %>% unique
allPwy=id_pwy$Pwy[!is.na(id_pwy$Pwy)] %>% unique
PwyToBeRemoved=allPwy[!(allPwy %in% PwyLeft)] #41 pwy to be removed

id_pwyLeft=id_pwy
id_pwyLeft$Pwy[id_pwyLeft$Pwy %in% PwyToBeRemoved]=NA #Note: NA %in% "A" gives FALSE
#------------------------------------------------------------------------------------


dat=All_Data_NAimputed
attribute_list=attr_list(id_pwyLeft$ids,id_pwyLeft$Pwy)

start.time = Sys.time()
corrVSAnnot=dist_TF_cumsum_matirxOperation(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
end.time = Sys.time()
end.time - start.time #Time difference of 39.01773 secs


sigPwy_confusionMatrix_metrics=confusionMatrix_metrics(corrVSAnnot$cumsum)


#================================================================================================



dist_list=list("Complete Pwy - PCC"=Original_pcc_confusionMatrix_metrics,
               "Significant Pwy - PCC"=sigPwy_confusionMatrix_metrics)
                     

##Set no. of samples to be plotted
samples=1:3000
##samples=seq(from=0,to=7914231,by=1989); samples[1]=1

p1=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                metric="sensitivity",y="Sensitivity") 
p2=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                metric="specificity",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                metric="precision",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                metric="accuracy",y="accuracy")


require(gridExtra)
p_all=grid.arrange(p1,p2,p3,p4, ncol=2) #When viewed in small window the dashed lines don't show properly




#Plot the above result using a fixed pcc. Plot: Coverage VS. accuracy + Coverage VS. sensitivity


###pick strain pairs that have PCC>=0.6 (1-|PCC|<=0.4) and are co-annotated by the same pathway(s)
cutoff=0.4 
dat=corrVSAnnot

strainPairs=dat[(dat$Distance<cutoff|dat$Distance==cutoff) & dat$sameORnot==1,c("Object 1","Object 2")]

allCoAnnotatedStrainPairs=dat[dat$sameORnot==1,c("Object 1","Object 2")]

capturedStrains_unique=c(strainPairs[,1],strainPairs[,2]) %>% unique

allCoAnnotatedStrains_unique=c(allCoAnnotatedStrainPairs[,1],allCoAnnotatedStrainPairs[,2]) %>% unique

###Retrieve the recovered pathway by the pcc cutoff
attribute_list=attr_list(id_pwyLeft$ids,id_pwyLeft$Pwy)

recoveredPwy=apply(strainPairs,1,FUN=function(strain1_strain2){
  strain1=as.character(strain1_strain2[1])
  strain2=as.character(strain1_strain2[2])
  
  intersect(attribute_list[[strain1]],attribute_list[[strain2]])
})

recoveredPwy_unique=unlist(recoveredPwy) %>% unique
str(recoveredPwy_unique)


#Get all possible IDs from the recovered pathway
allIDinRecoveredPwy=id_pwyLeft$ids[!is.na(id_pwyLeft[[2]]) & (id_pwyLeft[[2]] %in% recoveredPwy_unique)] %>% unique
str(allIDinRecoveredPwy)



#coverage = (IDs that are co-annotated in pathways by the cutoff) / (all possible IDs from the pathway co-annotations)
length(capturedStrains_unique)/length(allCoAnnotatedStrains_unique) #16.19%

#accuracy
index=sum((dat$Distance<cutoff|dat$Distance==cutoff))
confusionMatrix_metrics(corrVSAnnot$cumsum)$accuracy[index] #99.9% (I think it's so high because there are so many True negatives)

#sensitivity
confusionMatrix_metrics(corrVSAnnot$cumsum)$sensitivity[index]

#What's the percentage of pathways that are captured by such pcc cutoff?
length(recoveredPwy_unique)/length(id_pwyLeft$Pwy[!is.na(id_pwyLeft$Pwy)] %>% unique) #22.07%
  
  
#plot the result

#What would be the choose of parameters? => I think it can be different pcc cutoff
#What should be the variables?  Different combination of annotation sets  -> Type of data (There has to be different distance here) -> Distances

#What's the expectation?

  