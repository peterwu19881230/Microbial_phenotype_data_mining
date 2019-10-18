#Plot the above result using a fixed pcc. Plot: Coverage VS. accuracy + Coverage VS. sensitivity


#Calculate confusion matrix and other metrices using the annotations that have at least 1 sig phenotype
#================================================================================================


#source the function to get the confusion matrix and related stuff
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste(currentDir,"confusionMatrix_metrics.R",sep="/"))


#This part removes the annotations where non of the annotated genes have significant phenotypes
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste(currentDir,"removeNonSigAnnot.R",sep="/"))


corrVSannotANDmetrics=function(attribute_list){
  
  corrVSannot=dist_TF_cumsum_matirxOperation(data=All_Data_NAimputed,
                                             attribute_list=attr_list(attribute_list[[1]],attribute_list[[2]]), #[[1]] is the id, [[2]] is the annotation
                                             dist_metric=pcc_dist)
  
  metrics_df=confusionMatrix_metrics(corrVSannot$cumsum)
  
  return(list(corrVSannot=corrVSannot, metrics_df=metrics_df))
}



#start.time = Sys.time()

#pwy=corrVSannotANDmetrics(id_pwy)
#pcomplex=corrVSannotANDmetrics(id_pcomplex)
#pwyANDpcomplex=corrVSannotANDmetrics(id_pwyANDpcomplex)
#allAnnotation=corrVSannotANDmetrics(id_allAnnot)

#end.time = Sys.time()
#end.time - start.time  #Time difference of 13.08516 mins
#save(pwy,pcomplex,pwyANDpcomplex,allAnnotation,file="Data/sigAnnot_corrVSannot.RData")
load("Data/sigAnnot_corrVSannot.RData") #this takes a while because the files is 1.77G




#================================================================================================


###pick strain pairs that have PCC>=0.6 (1-|PCC|<=0.4) and are co-annotated by the same annotations(s)

metricByCutoff=function(cutoff=0.4,id_annot,corrVSAnnot_dat,metrics_df){
  
  dat=corrVSAnnot_dat
  attribute_list=attr_list(id_annot[[1]],id_annot[[2]])
  
  
  
  strainPairs=dat[(dat$Distance<cutoff|dat$Distance==cutoff) & dat$sameORnot==1,c("Object 1","Object 2")]
  
  allCoAnnotatedStrainPairs=dat[dat$sameORnot==1,c("Object 1","Object 2")]
  
  capturedStrains_unique=c(strainPairs[,1],strainPairs[,2]) %>% unique
  
  allCoAnnotatedStrains_unique=c(allCoAnnotatedStrainPairs[,1],allCoAnnotatedStrainPairs[,2]) %>% unique
  
  
  recoveredAnnot=apply(strainPairs,1,FUN=function(strain1_strain2){
    strain1=as.character(strain1_strain2[1])
    strain2=as.character(strain1_strain2[2])
    
    intersect(attribute_list[[strain1]],attribute_list[[strain2]])
  })
  recoveredAnnot_unique=unlist(recoveredAnnot) %>% unique
  
  
  recoveredAnnot_all=apply(allCoAnnotatedStrainPairs,1,FUN=function(strain1_strain2){
    strain1=as.character(strain1_strain2[1])
    strain2=as.character(strain1_strain2[2])
    
    intersect(attribute_list[[strain1]],attribute_list[[strain2]])
  })
  recoveredAnnot_all_unique=unlist(recoveredAnnot_all) %>% unique
  
  
  #Get all possible IDs from the recovered annotations
  allIDinRecoveredAnnot=id_annot$ids[!is.na(id_annot[[2]]) & (id_annot[[2]] %in% recoveredAnnot_unique)] %>% unique
  
  #coverage = (IDs that are co-annotated in annotations by the cutoff) / (all possible IDs from the co-annotations)
  coverage=length(capturedStrains_unique)/length(allCoAnnotatedStrains_unique) 
  
  #accuracy
  index=sum(dat$Distance<cutoff|dat$Distance==cutoff)
  accuracy=metrics_df$accuracy[index] #I think it would be so high because there are so many True negatives
  
  #sensitivity
  sensitivity=metrics_df$sensitivity[index]
  
  #precision
  precision=metrics_df$precision[index]
  
  #What's the percentage of unique annotations that are captured by such pcc cutoff?
  percentageUniqueAnnot=length((recoveredAnnot_unique))/length(recoveredAnnot_all_unique)
  
  
  return(c(coverage=coverage, accuracy=accuracy, sensitivity=sensitivity, precision=precision, percentageUniqueAnnot=percentageUniqueAnnot))
}


cutoffs=1-c(0.9,0.8,0.7,0.6,0.1)


pwyResult=metricByCutoff(cutoff=cutoffs[1],id_annot=id_pwy,corrVSAnnot_dat=pwy[[1]],metrics_df=pwy[[2]])
pwyResult=rbind(pwyResult,metricByCutoff(cutoff=cutoffs[2],id_annot=id_pwy,corrVSAnnot_dat=pwy[[1]],metrics_df=pwy[[2]]))
pwyResult=rbind(pwyResult,metricByCutoff(cutoff=cutoffs[3],id_annot=id_pwy,corrVSAnnot_dat=pwy[[1]],metrics_df=pwy[[2]]))
pwyResult=rbind(pwyResult,metricByCutoff(cutoff=cutoffs[4],id_annot=id_pwy,corrVSAnnot_dat=pwy[[1]],metrics_df=pwy[[2]])) #this result is consistent with the result in confusionMatrix_PwyWithSigFitness.R
pwyResult=rbind(pwyResult,metricByCutoff(cutoff=cutoffs[5],id_annot=id_pwy,corrVSAnnot_dat=pwy[[1]],metrics_df=pwy[[2]]))
pwyResult=data.frame(pwyResult,group="Pathway")


pcomplexResult=metricByCutoff(cutoff=cutoffs[1],id_annot=id_pcomplex,corrVSAnnot_dat=pcomplex[[1]],metrics_df=pcomplex[[2]]) 
pcomplexResult=rbind(pcomplexResult,metricByCutoff(cutoff=cutoffs[2],id_annot=id_pcomplex,corrVSAnnot_dat=pcomplex[[1]],metrics_df=pcomplex[[2]])) 
pcomplexResult=rbind(pcomplexResult,metricByCutoff(cutoff=cutoffs[3],id_annot=id_pcomplex,corrVSAnnot_dat=pcomplex[[1]],metrics_df=pcomplex[[2]])) 
pcomplexResult=rbind(pcomplexResult,metricByCutoff(cutoff=cutoffs[4],id_annot=id_pcomplex,corrVSAnnot_dat=pcomplex[[1]],metrics_df=pcomplex[[2]])) 
pcomplexResult=rbind(pcomplexResult,metricByCutoff(cutoff=cutoffs[5],id_annot=id_pcomplex,corrVSAnnot_dat=pcomplex[[1]],metrics_df=pcomplex[[2]])) 
pcomplexResult=data.frame(pcomplexResult,group="Protein complex")


pwyANDpcomplexResult=metricByCutoff(cutoff=cutoffs[1],id_annot=id_pwyANDpcomplex,corrVSAnnot_dat=pwyANDpcomplex[[1]],metrics_df=pwyANDpcomplex[[2]]) 
pwyANDpcomplexResult=rbind(pwyANDpcomplexResult,metricByCutoff(cutoff=cutoffs[2],id_annot=id_pwyANDpcomplex,corrVSAnnot_dat=pwyANDpcomplex[[1]],metrics_df=pwyANDpcomplex[[2]])) 
pwyANDpcomplexResult=rbind(pwyANDpcomplexResult,metricByCutoff(cutoff=cutoffs[3],id_annot=id_pwyANDpcomplex,corrVSAnnot_dat=pwyANDpcomplex[[1]],metrics_df=pwyANDpcomplex[[2]])) 
pwyANDpcomplexResult=rbind(pwyANDpcomplexResult,metricByCutoff(cutoff=cutoffs[4],id_annot=id_pwyANDpcomplex,corrVSAnnot_dat=pwyANDpcomplex[[1]],metrics_df=pwyANDpcomplex[[2]])) 
pwyANDpcomplexResult=rbind(pwyANDpcomplexResult,metricByCutoff(cutoff=cutoffs[5],id_annot=id_pwyANDpcomplex,corrVSAnnot_dat=pwyANDpcomplex[[1]],metrics_df=pwyANDpcomplex[[2]])) 
pwyANDpcomplexResult=data.frame(pwyANDpcomplexResult,group="Pathway & Protein complex")


allAnnotResult=metricByCutoff(cutoff=cutoffs[1],id_annot=id_allAnnot,corrVSAnnot_dat=allAnnotation[[1]],metrics_df=allAnnotation[[2]]) 
allAnnotResult=rbind(allAnnotResult,metricByCutoff(cutoff=cutoffs[2],id_annot=id_allAnnot,corrVSAnnot_dat=allAnnotation[[1]],metrics_df=allAnnotation[[2]])) 
allAnnotResult=rbind(allAnnotResult,metricByCutoff(cutoff=cutoffs[3],id_annot=id_allAnnot,corrVSAnnot_dat=allAnnotation[[1]],metrics_df=allAnnotation[[2]])) 
allAnnotResult=rbind(allAnnotResult,metricByCutoff(cutoff=cutoffs[4],id_annot=id_allAnnot,corrVSAnnot_dat=allAnnotation[[1]],metrics_df=allAnnotation[[2]])) 
allAnnotResult=rbind(allAnnotResult,metricByCutoff(cutoff=cutoffs[5],id_annot=id_allAnnot,corrVSAnnot_dat=allAnnotation[[1]],metrics_df=allAnnotation[[2]])) 
allAnnotResult=data.frame(allAnnotResult,group="All annotation sets")



#plot the result

#What would be the choose of parameters? => I think it can be different pcc cutoff
#What should be the variables?  Different combination of annotation sets  -> Type of data (There has to be different distance here) -> Distances

result=rbind(pwyResult,pcomplexResult,pwyANDpcomplexResult,allAnnotResult)
result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend


#sensitivity VS. precision (axes can be changed to have different metrics on x and y)
#(!)The result of removing annotations that don't have significant phenotypes look almost the same as if I don't remove them => Removing or not doesn't matter

Cols=c("#56B4E9","#F3518A","#BE70EA","#09D38A") 

p=ggplot(data=result,aes(sensitivity,precision,color=group))+

  geom_line()+
  geom_point(size=3)+
  geom_text(label=as.character(1-cutoffs) %>% rep(4),
            vjust=rep(-1,4*5),size=5,show.legend = FALSE)+ #show.legend = FALSE removes the mysterious "a" in the legend (https://stackoverflow.com/questions/18337653/remove-a-from-legend-when-using-aesthetics-and-geom-text)
  theme(plot.title= element_text(size = 20),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=20))+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  scale_color_manual(values=Cols)  
  

#currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
#pdf(file=paste(currentDir,"sensitivity_precision_sigAnnot.pdf",sep="/"),width=3*4,height=2*4)
#p
#dev.off()





