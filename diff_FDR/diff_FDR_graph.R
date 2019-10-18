#(!) The code here should be executable but computing MI takes lots of time, as apposed to pcc. So I will recode.

#Generate corr VS annot graph based on these parameters:
#1. Ternary data using different FDRs
#2. different annotations
#3. differnt metrics calculated from confusion matrix


#specify parameters
#====================================================================================================
#1. load various data using different FDR
load("Data/ternary_data_various_FDR.RData") #This gives: 
##     Ternary_Data_324cutff_NAremoved_1_FDR,
##     Ternary_Data_324cutff_NAremoved_10_FDR,
##     Ternary_Data_324cutff_NAremoved_15_FDR,
##     Ternary_Data_324cutff_NAremoved_20_FDR,

##Ternary_Data_324cutff_NAremoved is already sourced from Nichols_preload.R. This one uses 5% FDR

#2.
annots=c("Pwy","pcomplex","operon","regulator","kegg_modules") #("any" will be constructed later)

#3.
metrics=c("precision","sensitivity")
#====================================================================================================


#prep data
#====================================================================================================
my_mutual_info_dist=function(dat){
  mi=infotheo::natstobits(infotheo::mutinformation(as.data.frame(t(dat))))
  mi_distObj=as.dist(mi)
  mi_distance=1-mi_distObj
  return(mi_distance)
}


#Parallel computing material
################
library(foreach)
library(doParallel)
library(parallel)
library(doSNOW) #for windows I have to use this

numCores=detectCores()-1
cl=makeCluster(numCores)
registerDoSNOW(cl) 

################


start.time=Sys.time()



ternary_list=list(Ternary_Data_324cutff_NAremoved_1_FDR,
                  Ternary_Data_324cutff_NAremoved,
                  Ternary_Data_324cutff_NAremoved_10_FDR,
                  Ternary_Data_324cutff_NAremoved_15_FDR,
                  Ternary_Data_324cutff_NAremoved_20_FDR)

result=foreach(dat=ternary_list) %dopar%{
  my_mutual_info_dist(dat)
}


stopCluster(cl)
end.time=Sys.time()
end.time-start.time #Time difference of 1.815686 hours on my PC


pairwise_mi_diffCutoff=result

#save(pairwise_mi_diffCutoff,file="Data/pairwise_mi_diffCutoff.RData")


pairwise_mi_diffCutoff_dfs=lapply(pairwise_mi_diffCutoff,FUN=melt_dist)


pairwise_mi_diffCutoff_completeDF=cbind(pairwise_mi_diffCutoff_dfs[[1]],
                                        pairwise_mi_diffCutoff_dfs[[2]][,3],
                                        pairwise_mi_diffCutoff_dfs[[3]][,3],
                                        pairwise_mi_diffCutoff_dfs[[4]][,3],
                                        pairwise_mi_diffCutoff_dfs[[5]][,3])

names(pairwise_mi_diffCutoff_completeDF)=c("strain1","strain2","1_FDR","5_FDR","10_FDR","15_FDR","20_FDR")
pairwise_mi_diffCutoff_completeDF$strain1=as.numeric(pairwise_mi_diffCutoff_completeDF$strain1)
pairwise_mi_diffCutoff_completeDF$strain2=as.numeric(pairwise_mi_diffCutoff_completeDF$strain2)


mi_diffCutoff_completePairwiseTable=merge(pairwise_mi_diffCutoff_completeDF,
                                          strain1strain2_allAnnotations_allDistances[,c(c("strain1","strain2"),annots)],
                                          by=c("strain1","strain2"))

mi_diffCutoff_completePairwiseTable[,"any"]=ifelse(rowSums(mi_diffCutoff_completePairwiseTable[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>0,1,0)


#save(mi_diffCutoff_completePairwiseTable,file="Data/mi_diffCutoff_completePairwiseTable.RData")

load("Data/mi_diffCutoff_completePairwiseTable.RData")

confusionMatrix_metrics_list=list()
for(FDR in c("1_FDR","5_FDR","10_FDR","15_FDR","20_FDR")){
  
  for(annot in c("Pwy","pcomplex","operon","regulator","kegg_modules","any")){
    cumsum_=mi_diffCutoff_completePairwiseTable[order(mi_diffCutoff_completePairwiseTable[,FDR]),annot] %>% cumsum
    confusionMatrix_metrics_list[[FDR]][[annot]]=confusionMatrix_metrics(cumsum_)
  }
  
}


confusionMatrix_metrics_list=purrr::transpose(confusionMatrix_metrics_list) #turn the list inside-out


#====================================================================================================


#get graphs
#====================================================================================================
plot_metrics=function(samples,alpha,size,annot_list,metric,x_lab="high similarity -- ranked pairs -- low similarity",y,legend_title="Annotation set(s)"){
  
  random_metric=paste("random",metric,sep="_")
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(no=samples,Metric=annot_list[[1]][[metric]][samples]),aes(no,Metric,color=names(annot_list[1])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=annot_list[[2]][[metric]][samples]),aes(no,Metric,color=names(annot_list[2])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=annot_list[[3]][[metric]][samples]),aes(no,Metric,color=names(annot_list[3])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=annot_list[[4]][[metric]][samples]),aes(no,Metric,color=names(annot_list[4])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=annot_list[[5]][[metric]][samples]),aes(no,Metric,color=names(annot_list[5])),size=size)+
    
    
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=annot_list[[1]][[random_metric]][samples]),aes(no,Metric,color=names(annot_list[1])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=annot_list[[2]][[random_metric]][samples]),aes(no,Metric,color=names(annot_list[2])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=annot_list[[3]][[random_metric]][samples]),aes(no,Metric,color=names(annot_list[3])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=annot_list[[4]][[random_metric]][samples]),aes(no,Metric,color=names(annot_list[4])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=annot_list[[5]][[random_metric]][samples]),aes(no,Metric,color=names(annot_list[5])),size=size,alpha=alpha,linetype = "dashed")+
    
    
    
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)
}

samples_all=seq(from=0,to=7914231,by=1989); samples_all[1]=1

samples_list=list(part=1:4000,all=samples_all) 
annots=c("Pwy","pcomplex","operon","regulator","kegg_modules","any")
metrics=c("sensitivity","precision")


dir_of_workingScript=dirname(rstudioapi::getActiveDocumentContext()$path)

start.time=Sys.time()

for(annot in annots){
  for(metric in metrics){
    for(samples in samples_list){
      p=plot_metrics(samples=samples,alpha=0.5,size=1,annot_list=confusionMatrix_metrics_list[[annot]],metric=metric,y=metric)
      
      ggsave(file=paste(dir_of_workingScript,"/","mi","_",annot,"_diff_FDR_",metric,"_",max(samples),".png",sep=""),plot=p,device="png",width=3*4,height=2*4)
      }
    
  }
  
}

end.time=Sys.time()
end.time-start.time #Time difference of 42.13701 secs
#====================================================================================================



