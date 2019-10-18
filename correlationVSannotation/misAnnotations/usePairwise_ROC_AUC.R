#use pairwise data to generate mis-annotations and plot ROC, calculate AUC

df=strain1strain2_allAnnotations_allDistances[order(strain1strain2_allAnnotations_allDistances$pcc),c("Pwy","pcomplex","operon","regulator","kegg_modules")]

coAnnotation=ifelse(rowSums(df)==5,1,0)
sum(coAnnotation)

#Function to generate a df that has original annotations + different no. of annotations left out
remove_coAnnotation=function(coAnnotation,fraction=c(0.05,0.25,0.5,0.6,0.7,0.8,0.9)){
  
  leng=length(coAnnotation)
  tot=sum(coAnnotation)
  n=round(tot*fraction)
  IndexForCoAnnotation=which(coAnnotation==1)
  
  
  new_coAnnotation=coAnnotation
  for(i in seq(n)){
    
    indexToBeRemoved=sample(IndexForCoAnnotation,n[i])
    new=coAnnotation
    new[indexToBeRemoved]=0
    
    print(sum(coAnnotation!=new)) #confirm the no. of mis-annotations
    
    new_coAnnotation=cbind(new_coAnnotation,new)
    
  }
  
  
  colnames(new_coAnnotation)=c("original",fraction) 
  
  return(new_coAnnotation)
}

#1. Function to generate a df that has original annotations + different no. of shuffled annotations
#Note: the no. of misannotation doubles the provided "fraction" because 1 shuffling gives 2 annotations (1 missed + 1 added)
shuffle_coAnnotation=function(coAnnotation,fraction=c(0.05,0.25,0.5,0.6,0.7,0.8,0.9)){
  
  leng=length(coAnnotation)
  tot=sum(coAnnotation)
  n=round(tot*fraction)
  IndexForCoAnnotation=which(coAnnotation==1)
  
  
  new_coAnnotation=coAnnotation
  for(i in seq(n)){
    
    indexToBeRemoved=sample(IndexForCoAnnotation,n[i])
    indexToBePut=sample((1:leng)[-IndexForCoAnnotation],n[i])
    
    new=coAnnotation
    new[indexToBeRemoved]=0
    new[indexToBePut]=1
    
    
    print(sum(coAnnotation!=new)) #confirm the no. of mis-annotations
    
    new_coAnnotation=cbind(new_coAnnotation,new)
    
  }
  
  colnames(new_coAnnotation)=c("original",fraction) 
  
  return(new_coAnnotation)
}


removed_cumsum=remove_coAnnotation(coAnnotation) %>% apply(2,cumsum)
shuffled_cumsum=shuffle_coAnnotation(coAnnotation) %>% apply(2,cumsum)


removed_confuion_metrics=apply(removed_cumsum,2,confusionMatrix_metrics)
shuffled_confuion_metrics=apply(shuffled_cumsum,2,confusionMatrix_metrics)




#plot ROC
plot_ROC_fromCoAnnotation=function(samples,alpha,size,conf_metrics_list,legend_title="Mis-annotation set(s)"){
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  p=ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                   legend.text=element_text(size=15),
                                   axis.text.x = element_text(size=15),
                                   axis.text.y=element_text(size=15),
                                   axis.title=element_text(size=20))+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[1]]$specificity[samples],sensitivity=conf_metrics_list[[1]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[1])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[2]]$specificity[samples],sensitivity=conf_metrics_list[[2]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[2])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[3]]$specificity[samples],sensitivity=conf_metrics_list[[3]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[3])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[4]]$specificity[samples],sensitivity=conf_metrics_list[[4]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[4])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[5]]$specificity[samples],sensitivity=conf_metrics_list[[5]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[5])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[6]]$specificity[samples],sensitivity=conf_metrics_list[[6]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[6])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[7]]$specificity[samples],sensitivity=conf_metrics_list[[7]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[7])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[8]]$specificity[samples],sensitivity=conf_metrics_list[[8]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[8])),size=size)+
    
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagonal line
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)
  
  return(p)
}


samples=seq(from=0,to=7914231,by=1989); samples[1]=1
roc_remove=plot_ROC_fromCoAnnotation(samples=samples,alpha=0.8,size=1,removed_confuion_metrics)
roc_shuffle=plot_ROC_fromCoAnnotation(samples=samples,alpha=0.8,size=1,shuffled_confuion_metrics)


roc_remove
roc_shuffle 
 

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"roc_remove.png",sep="/"),plot=roc_remove,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension



dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"roc_shuffle.png",sep="/"),plot=roc_shuffle,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension











