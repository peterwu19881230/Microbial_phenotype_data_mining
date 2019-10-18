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




#Throw out conditions that don't have many significant phenotypes for genes involved in EcoCyc pathways
#============================================================================================================================================
index=id_allAttributes$ids[!is.na(id_allAttributes$Pwy)] %>% unique
length(index) 
##note: index is a character vector, but it can be used to subset the following data which has the character rownames

datInPwy=Ternary_Data_324cutff_NAremoved[index,]
dim(datInPwy)

noOfSigPhentypes=apply(Ternary_Data_324cutff_NAremoved,2,FUN=function(col) sum(col!=0))

fitness=as.numeric(datInPwy)
success=sum(fitness!=0) #6318

#perform over- or under-representation (Ref: http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html)
#Should I use over- or under- ?

#over-representation
hitInSample=noOfSigPhentypes #noOfSigPhentypes is defined above
hitInPop=success
failInPop=dim(datInPwy)[1]*dim(datInPwy)[2]-success
sampleSize=dim(datInPwy)[1]


pVals=phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
sigPVals=pVals[pVals<=0.05]
length(sigPVals) #218 conditions to be used
#============================================================================================================================================




#Calculate confusion matrix and other metrices after using new conditions
#================================================================================================
condionIndex=(pVals<=0.05) 

dat=All_Data_NAimputed[,condionIndex]
id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique()
attribute_list=attr_list(id_pwy$ids,id_pwy$Pwy)

start.time = Sys.time()
corrVSAnnot=dist_TF_cumsum_matirxOperation(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
end.time = Sys.time()
end.time - start.time #Time difference of 37.16093 secs


pcc_confusionMatrix_metrics=confusionMatrix_metrics(corrVSAnnot$cumsum)
#================================================================================================



dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
               "Data of enriched conditions - PCC"=pcc_confusionMatrix_metrics)
                     

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




#Plot the above result using a fixed pcc. Plot: Coverage VS. accuracy


###pick strain pairs that have PCC>=0.6 (1-|PCC|<=0.4) and are co-annotated by the same pathway(s)
cutoff=0.4 
dat=strain1strain2_allAnnotations_allDistances
strainPairs=dat[(dat$pcc<cutoff|dat$pcc==cutoff) & dat$Pwy==1,c("strain1","strain2")]

capturedStrains_unique=c(strainPairs[,1],strainPairs[,2]) %>% unique


###Retrieve the recovered pathway by the pcc cutoff
id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique()
attribute_list=attr_list(id_pwy$ids,id_pwy$Pwy)

recoveredPwy=apply(strainPairs,1,FUN=function(strain1_strain2){
  strain1=as.character(strain1_strain2[1])
  strain2=as.character(strain1_strain2[2])
  
  intersect(attribute_list[[strain1]],attribute_list[[strain2]])
})

recoveredPwy_unique=unlist(recoveredPwy) %>% unique
str(recoveredPwy_unique)


#Get all possible IDs from the recovered pathway
allIDinRecoveredPwy=id_allAttributes$ids[!is.na(id_allAttributes$Pwy) & id_allAttributes$Pwy %in% recoveredPwy_unique] %>% unique
str(allIDinRecoveredPwy)

#coverage = (IDs that are co-annotated in pathways by the cutoff) / (all possible IDs from the recovered pathway)
length(capturedStrains_unique)/length(allIDinRecoveredPwy) #36.2%

#accuracy
index=sum((dat$pcc<cutoff|dat$pcc==cutoff))
Original_pcc_confusionMatrix_metrics$accuracy[index] #99.9% (I think it's so high because there are so many True negatives)


#What's the percentage of pathways that are captured by such pcc cutoff?
length((recoveredPwy_unique))/length(id_allAttributes$Pwy %>% unique) #20.1%
  
  
#plot the result

#What would be the choices of parameters?
#What should be the variables?  Type of -> Type of data (There has to be different distance here) -> Distances

#What's the expectation?

  