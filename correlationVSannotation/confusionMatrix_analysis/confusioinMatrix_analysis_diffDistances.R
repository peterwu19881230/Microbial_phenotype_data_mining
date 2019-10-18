#x: ranked pairs , y: Sensitivity --- different combinatio of annotation sets on the same graph
#(From confusionMatrix_analysis.R I noticed that sensitivuty is probably the only thing that is worth plotting)



#4 metrics defined in Wikipedia:
##sensitivity (recall) = TP / (TP +FN)
##specificity = TN / (TN+FP)
##precision = TP/(TP+FP)
##accuracy = (TP + TN)/total


##suppose the cutoff pcc=rho

##TP=no. of pairs that: (a) have pcc >= rho (b) are co-annotated
##TN=no. of pairs that: (a) have pcc < rho (b) are not co-annotated
##FP=no. of pairs that: (a) have pcc >= rho (b) are not co-annotated
##FN=no. of pairs that: (a) have pcc < rho (b) are co-annotated





#Function to get the graph.
#========================================================================================================================
plot_metrics=function(samples,alpha,size,dist_list,metric,
                      index_collapsed_MHD3=NULL, #These are to give indices where rankings make real sense
                      x_lab="Low distance -- ranked pairs -- High distance",y,legend_title="Distance"){
  
  random_metric=paste("random",metric,sep="_")
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
    geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],Metric=dist_list[[1]][[metric]][index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,Metric,color=names(dist_list[1])),size=size*3,alpha=alpha,shape=16)+
     geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],Metric=dist_list[[1]][[metric]][index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,Metric,color=names(dist_list[1])),size=size,alpha=alpha)+
    
    geom_line(data=data.frame(no=samples,Metric=dist_list[[2]][[metric]][samples]),aes(no,Metric,color=names(dist_list[2])),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[3]][[metric]][samples]),aes(no,Metric,color=names(dist_list[3])),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[4]][[metric]][samples]),aes(no,Metric,color=names(dist_list[4])),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[5]][[metric]][samples]),aes(no,Metric,color=names(dist_list[5])),size=size,alpha=alpha)+
    
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=dist_list[[1]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[1])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[2]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[2])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[3]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[3])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[4]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[4])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[5]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[5])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    
    
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)
}
#========================================================================================================================



#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc","mhd3_collapsedCond","mi_ternary","mi","spearman")]
TF=( rowSums(df[,-c(6,7,8,9,10)])==5 )
sum(TF) #106


pcc=cumsum( TF[order(df$pcc)])
pcc_confusionMatrix_metrics=confusionMatrix_metrics(pcc)

spearman=cumsum( TF[order(df$spearman)])
spearman_confusionMatrix_metrics=confusionMatrix_metrics(spearman)

mi=cumsum( TF[order(df$mi)])
mi_confusionMatrix_metrics=confusionMatrix_metrics(mi)


TernaryCondCollapsed_mhd3=cumsum( TF[order(df$mhd3_collapsedCond)])
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt
TernaryCondCollapsed_mhd3_confusionMatrix_metrics=confusionMatrix_metrics(TernaryCondCollapsed_mhd3)


TernaryCondCollapsed_mi=cumsum( TF[order(df$mi_ternary)])
TernaryCondCollapsed_mi_confusionMatrix_metrics=confusionMatrix_metrics(TernaryCondCollapsed_mi)




dist_list=list(
               "MHD"=TernaryCondCollapsed_mhd3_confusionMatrix_metrics,
               "pcc"=pcc_confusionMatrix_metrics,
               "spearman"=spearman_confusionMatrix_metrics,
               "MI"=mi_confusionMatrix_metrics,
               "Ternary MI"=TernaryCondCollapsed_mi_confusionMatrix_metrics
                   )





##Set no. of samples to be plotted
##samples=1:300
samples=1:5000
##samples=seq(from=0,to=7914231,by=1989); samples[1]=1


index_collapsed_MHD3=df$mhd3_collapsedCond %>% table %>% cumsum

p1=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                index_collapsed_MHD3,
                metric="sensitivity",y="Sensitivity") 
p2=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                index_collapsed_MHD3,
                metric="specificity",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                index_collapsed_MHD3,
                metric="precision",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                index_collapsed_MHD3,
                metric="accuracy",y="accuracy")


require(gridExtra)
p_all=grid.arrange(p1,p2,p3,p4, ncol=2) #When viewed in small window the dashed lines don't show properly

#dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
#ggsave(file=paste(dir_of_workingScript,"metrics_diffDistance.pdf",sep="/"),plot=p_all,device="pdf",width=3*6,height=2*6)

##When samples=1:300 or =1:5000 it cuts the edges so I changed the width and hieight a little bit
#ggsave(file=paste(dir_of_workingScript,"metrics_diffDistance_first5000.pdf",sep="/"),plot=p_all,device="pdf",width=3*7,height=2*7)
#ggsave(file=paste(dir_of_workingScript,"metrics_diffDistance_first300.pdf",sep="/"),plot=p_all,device="pdf",width=3*7,height=2*7)






#sensitivity VS. precision (axes can be changed to have different metrics on x and y)

#Complete graph sampling from all 7914231 pairs
#============================================================================================================================================


metricByCutoff=function(metric_list,samples){
  
  result_list=list()
  
  for(i in seq(metric_list)){
    
    
    result_list[[i]]=data.frame(sensitivity=metric_list[[i]]$sensitivity[samples],
                                specificity=metric_list[[i]]$specificity[samples],
                                precision=metric_list[[i]]$precision[samples],
                                accuracy=metric_list[[i]]$accuracy[samples],
                                group=names(metric_list)[i])
    
  }
  
  result=Reduce(rbind,result_list) 
  result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend
  
  return(result)
}

samples=seq(from=0,to=7914231,by=1989)
samples=samples[-1] #throw out the first index because it would cause problem when doing log() on plotting
samples=c(2,seq(from=5,to=1989,by=5),samples) #I don't want to miss points using stringent cutoffs (ranking of of pairs <1989)

result=metricByCutoff(dist_list,samples)



Cols=c("black","#09D38A","#564c99","#832B05","#EC8D3A") 


p=ggplot(data=result,aes(sensitivity,precision,color=group))+
  theme_minimal()+
  geom_line()+
  #geom_point(size=3)+
  theme(plot.title= element_text(size = 20),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=20))+
  scale_y_continuous(limits = c(0,1),labels = function(x) paste0(round(x*100,3), "%"))+ 
  scale_x_continuous(limits = c(0,1),labels = function(x) paste0(round(x*100,3), "%"))+ 
  scale_color_manual(values=Cols)+
  labs(color="Type of data - Correlation")+
  coord_flip() #this flips x and y axes
#============================================================================================================================================

#A graph that only shows good cutoffs
#============================================================================================================================================

#Find the index where each point has the highest precision
optimal_indices=c()
for(i in seq(dist_list)){
  optimal_indices[i]=which(dist_list[[i]]$precision==max(dist_list[[i]]$precision))
}


metricByCutoff_2=function(metric_list,optimal_indices){
  
  
  result_list=list()
  
  for(i in seq(metric_list)){
    
    
    result_list[[i]]=data.frame(sensitivity=metric_list[[i]]$sensitivity[optimal_indices[[i]]],
                                specificity=metric_list[[i]]$specificity[optimal_indices[[i]]],
                                precision=metric_list[[i]]$precision[optimal_indices[[i]]],
                                accuracy=metric_list[[i]]$accuracy[optimal_indices[[i]]],
                                group=names(metric_list)[i])
    
  }
  
  result=Reduce(rbind,result_list) 
  result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend
  
  return(result)
}

result_2=metricByCutoff_2(dist_list,optimal_indices)


#Find the point for Ternary MI that has the highest precision when sensitivity > 0.125
precision=dist_list$`Ternary MI`$precision
sensitivity=dist_list$`Ternary MI`$sensitivity
Index=which(precision==max(precision[sensitivity>0.125])) #1680

result_2=rbind(result_2,data.frame(sensitivity=dist_list$`Ternary MI`$sensitivity[Index],
                                   specificity=dist_list$`Ternary MI`$specificity[Index],
                                   precision=dist_list$`Ternary MI`$precision[Index],
                                   accuracy=dist_list$`Ternary MI`$accuracy[Index],
                                   group="Ternary MI"))


#Find the point for spearman that has the highest precision when sensitivity > 0.2
precision=dist_list$`spearman`$precision
sensitivity=dist_list$`spearman`$sensitivity
Index=which(precision==max(precision[sensitivity>0.2])) #128

result_2=rbind(result_2,data.frame(sensitivity=dist_list$`spearman`$sensitivity[Index],
                                   specificity=dist_list$`spearman`$specificity[Index],
                                   precision=dist_list$`spearman`$precision[Index],
                                   accuracy=dist_list$`spearman`$accuracy[Index],
                                   group="spearman"))


result_2=result_2[-c(3,5),] #remove the non-optimal points


Cols=c("black","#09D38A","#564c99","#832B05","#EC8D3A") 


p=ggplot(data=result_2,aes(sensitivity,precision,color=group))+
  theme_minimal()+
  #geom_line()+
  geom_point(size=3)+
  theme(plot.title= element_text(size = 20),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=20))+
  scale_y_continuous(trans='log10',limits = c(0.001,1),labels = function(x) paste0(round(x*100,3), "%"))+ 
  scale_x_continuous(trans='log10',limits = c(0.001,1),labels = function(x) paste0(round(x*100,3), "%"))+ 
  scale_color_manual(values=Cols)+
  labs(color="Type of data - Correlation")+
  coord_flip() #this flips x and y axes


#The corresponding similarity
df=strain1strain2_allAnnotations_allDistances[,c("mhd3","pcc","spearman","mi","mi_ternary")]

new_optimal_indices=optimal_indices
new_optimal_indices[5]=1680 
new_optimal_indices[3]=128 

sim=c()
for(i in 1:dim(df)[2]){
  sim[i]=sort(df[,i])[new_optimal_indices[i]]
}

sim 
##MHD: 308.0000000   
##1 - |PCC|: 0.3277587      => |PCC| = 0.67
##1 - |Spearman|: 0.2774177 => |Spearman| = 0.72
##1 - MI: 0.4963942         => MI = 0.50   
##1- MI ternary: 0.8601969  => MI ternary = 0.14

#============================================================================================================================================


##For actual fig on the paper (1. not using MHD 2. Ternary MI is shown on another graph so I am not using it here)
#============================================================================================================================================
result_paper=result_2[!(result_2$group %in% c("MHD","Ternary MI")),]

Cols=c("#09D38A","#564c99","#832B05") 

p=ggplot(data=result_paper,aes(sensitivity,precision,color=group))+
  theme_minimal()+
  #geom_line()+
  geom_point(size=3)+
  theme(plot.title= element_text(size = 20),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=20))+
  scale_y_continuous(trans='log10',limits = c(0.01,1),labels = function(x) paste0(round(x*100,3), "%"))+ 
  scale_x_continuous(trans='log10',limits = c(0.01,1),labels = function(x) paste0(round(x*100,3), "%"))+ 
  scale_color_manual(values=Cols)+
  labs(color="Type of data - Correlation")+
  coord_flip() #this flips x and y axes

#============================================================================================================================================



currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
pdf(file=paste(currentDir,"sensitivity_precision_Dist.pdf",sep="/"),width=3*4,height=2*4)
p
dev.off()


# ROC-AUC
#============================================================================================================================================
plot_ROC_original=function(samples,alpha,size,dist_list,legend_title="Similarity"){
  
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("black","#832B05","#09D38A","#564c99","#EC8D3A")
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(FPR=1-dist_list[[1]]$specificity[samples],sensitivity=dist_list[[1]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[1])),size=size)+
    geom_line(data=data.frame(FPR=1-dist_list[[2]]$specificity[samples],sensitivity=dist_list[[2]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[2])),size=size)+
    geom_line(data=data.frame(FPR=1-dist_list[[3]]$specificity[samples],sensitivity=dist_list[[3]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[3])),size=size)+
    geom_line(data=data.frame(FPR=1-dist_list[[4]]$specificity[samples],sensitivity=dist_list[[4]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[4])),size=size)+
    geom_line(data=data.frame(FPR=1-dist_list[[5]]$specificity[samples],sensitivity=dist_list[[5]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[5])),size=size)+
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagnal line
    #scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    #scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)+
    scale_color_manual(values=Cols) 
}

samples=1:7914231
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1
roc_diffDist_original=plot_ROC(samples=samples,alpha=0.5,size=1,dist_list)


##Thise only plots pcc, spearman, MI
plot_ROC_paper=function(samples,alpha,size,dist_list,legend_title="Similarity"){ 
  
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#832B05","#09D38A","#564c99")
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    #geom_line(data=data.frame(FPR=1-dist_list[[1]]$specificity[samples],sensitivity=dist_list[[1]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[1])),size=size)+
    geom_line(data=data.frame(FPR=1-dist_list[[2]]$specificity[samples],sensitivity=dist_list[[2]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[2])),size=size)+
    geom_line(data=data.frame(FPR=1-dist_list[[3]]$specificity[samples],sensitivity=dist_list[[3]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[3])),size=size)+
    geom_line(data=data.frame(FPR=1-dist_list[[4]]$specificity[samples],sensitivity=dist_list[[4]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[4])),size=size)+
    #geom_line(data=data.frame(FPR=1-dist_list[[5]]$specificity[samples],sensitivity=dist_list[[5]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(dist_list[5])),size=size)+
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagnal line
    #scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    #scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)+
    scale_color_manual(values=Cols) 
}

samples=1:7914231
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1
roc_diffDist=plot_ROC_paper(samples=samples,alpha=0.5,size=1,dist_list)


currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)

start.time=Sys.time()
ggsave(file=paste(currentDir,"roc_diffDist.png",sep="/"),plot=roc_diffDist,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension
end.time=Sys.time()
end.time-start.time #Time difference of 3.307932 mins


#Calculate AUC (Area under curve) using a package
library(pROC) #ref: https://stackoverflow.com/questions/4903092/calculate-auc-in-r

df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")]
response=( rowSums(df)==5 )


#pcc
predictor=strain1strain2_allAnnotations_allDistances$pcc

start.time = Sys.time()
pccAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time #Time difference of 22.78993 secs

pccAUC #Area under the curve: 0.9131


#spearman
predictor=strain1strain2_allAnnotations_allDistances$spearman

start.time = Sys.time()
spearmanAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 

spearmanAUC #Area under the curve: 0.8938


#MI
predictor=strain1strain2_allAnnotations_allDistances$mi

start.time = Sys.time()
miAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 

miAUC #Area under the curve: 0.8968

#ternary MI
predictor=strain1strain2_allAnnotations_allDistances$mi_ternary

start.time = Sys.time()
ter_miAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 

ter_miAUC #Area under the curve: 0.8432


#ternary mhd3
predictor=strain1strain2_allAnnotations_allDistances$mhd3

start.time = Sys.time()
mhd3AUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 


mhd3AUC #Area under the curve: 0.8546
#============================================================================================================================================




#PR curve (x: recall(sensitivity), y: precision  ) (ref: https://classeval.wordpress.com/introduction/introduction-to-the-precision-recall-plot/)

samples=1:7914231
#samples=seq(from=0,to=7914231,by=173); samples[1]=1 #7914231=3^2*13*17*23*173

plot_PRC=function(samples,size,dist_list,legend_title="Type of data - Correlation"){
  
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#832B05","#09D38A","#564c99")
  
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    #geom_line(data=data.frame(recall=1-dist_list[[1]]$sensitivity[samples],precision=dist_list[[1]]$precision[samples]),aes(recall,precision,color=names(dist_list[1])),size=size)+
    geom_line(data=data.frame(recall=1-dist_list[[2]]$sensitivity[samples],precision=dist_list[[2]]$precision[samples]),aes(recall,precision,color=names(dist_list[2])),size=size)+
    geom_line(data=data.frame(recall=1-dist_list[[3]]$sensitivity[samples],precision=dist_list[[3]]$precision[samples]),aes(recall,precision,color=names(dist_list[3])),size=size)+
    geom_line(data=data.frame(recall=1-dist_list[[4]]$sensitivity[samples],precision=dist_list[[4]]$precision[samples]),aes(recall,precision,color=names(dist_list[4])),size=size)+
    #geom_line(data=data.frame(recall=1-dist_list[[5]]$sensitivity[samples],precision=dist_list[[5]]$precision[samples]),aes(recall,precision,color=names(dist_list[5])),size=size)+
    
    #this is the negative control
    geom_abline(slope=0,intercept=(dist_list[[1]]$TP[1] + dist_list[[1]]$FN[1])/dim(dist_list[[1]])[1],linetype="dashed")+ #intercept= Positive/(Negative + Positive) = Positive/Total = (TP + FN)/total
    
    
    #scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    #scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="Recall (sensitivity)",y="Precision",aesthetic='custom text',color=legend_title)+
    #scale_x_continuous(expand = c(0,0),labels=comma)+
    scale_color_manual(values=Cols) 
}


prc_diffDistance=plot_PRC(samples,size=1,dist_list)

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)

start.time=Sys.time()
ggsave(file=paste(dir_of_workingScript,"prc_diffDistance.png",sep="/"),plot=prc_diffDistance,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension
end.time=Sys.time()
end.time-start.time #Time difference of 4.063097 mins



library(PRROC)

df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")]
response=( rowSums(df)==5 )

#pcc
predictor=1-strain1strain2_allAnnotations_allDistances$pcc

pr.curve(scores.class0 = predictor,weights.class0=as.numeric(response))
##Area under curve (Integral): 0.01240156 
##Area under curve (Davis & Goadrich):0.01240133 


#MI
predictor=1-strain1strain2_allAnnotations_allDistances$mi

pr.curve(scores.class0 = predictor,weights.class0=as.numeric(response))
##Area under curve (Integral): 0.05286958 
##Area under curve (Davis & Goadrich): 0.05286234 


#spearman
predictor=1-strain1strain2_allAnnotations_allDistances$spearman

pr.curve(scores.class0 = predictor,weights.class0=as.numeric(response))
##Area under curve (Integral): 0.07894991 
##Area under curve (Davis & Goadrich): 0.0789306


#random (Positive/Total) * 1 :
106/7914231 #1.34e-05




