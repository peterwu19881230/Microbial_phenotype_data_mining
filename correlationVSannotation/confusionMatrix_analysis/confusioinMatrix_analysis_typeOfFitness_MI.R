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
plot_metrics=function(samples,alpha,size,dist_list,metric,x_lab="Low distance -- ranked pairs -- High distance",y,legend_title="Type of data - Correlation"){
  
  #Cols=c("#09D38A","#EC8D3A","#832B05")
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#832B05","#09D38A","#EC8D3A") 
  
  random_metric=paste("random",metric,sep="_")
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(no=samples,Metric=dist_list[[1]][[metric]][samples]),aes(no,Metric,color=names(dist_list[1])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[2]][[metric]][samples]),aes(no,Metric,color=names(dist_list[2])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[3]][[metric]][samples]),aes(no,Metric,color=names(dist_list[3])),size=size)+
    
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=dist_list[[1]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[1])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[2]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[2])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=dist_list[[3]][[random_metric]][samples]),aes(no,Metric,color=names(dist_list[3])),size=size,alpha=alpha,linetype = "dashed")+
    
    
    
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)+
    scale_color_manual(values=rep(Cols,2)) 
}
#========================================================================================================================

inPrice=strain1strain2_allAnnotations_allDistances$inPriceOrNot

#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[inPrice,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc","mi_ternary","mi_ternary_collapsedCond")]
TF=( rowSums(df[,-c(6,7,8)])==5 )
sum(TF) #100



Original_pcc=cumsum( TF[order(df$pcc)])
Original_pcc_confusionMatrix_metrics=confusionMatrix_metrics(Original_pcc,total=6211050)

#start here. correct the following 2 lines so it uses intersection of all annotation sets

load("Data/Price_pairwisePCC_cumsum.RData")
Price_pcc_confusionMatrix_metrics=confusionMatrix_metrics(Price_pairwisePCC_cumsum$cumsum_all,total=6211050)


mi_ternary=cumsum( TF[order(df$mi_ternary)]) 
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt with
mi_ternary_confusionMatrix_metrics=confusionMatrix_metrics(mi_ternary,total=6211050)


mi_ternary_collapsedCond=cumsum( TF[order(df$mi_ternary_collapsedCond)])
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt with
mi_ternary_collapsedCond_confusionMatrix_metrics=confusionMatrix_metrics(mi_ternary_collapsedCond,total=6211050)


load("Data/NichDeut_mi_ternary_allAnnot.RData") 
NichDeut_mi_ternary_confusionMatrix_metrics=confusionMatrix_metrics(NichDeut_mi_ternary_allAnnot$cumsum,total=6211050)

load("Data/price_mi_ternary_allAnnot.RData")
Price_mi_ternary_confusionMatrix_metrics=confusionMatrix_metrics(price_mi_ternary_allAnnot$cumsum,total=6211050)

#dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
#                     "Ternary data - MI"=mi_ternary_confusionMatrix_metrics,
#                     "Ternary data with condition collapsed - MI"=mi_ternary_collapsedCond_confusionMatrix_metrics)

dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
               "Price original data - PCC"=Price_pcc_confusionMatrix_metrics,
               "Ternary data - MI"=mi_ternary_confusionMatrix_metrics,
               "Price Ternary data - MI"=Price_mi_ternary_confusionMatrix_metrics,
               "2 combined ternary data - MI"=NichDeut_mi_ternary_confusionMatrix_metrics)

##Set no. of samples to be plotted
##samples=1:300
##samples=1:4000
##samples=seq(from=0,to=7914231,by=1989); samples[1]=1
#samples=seq(from=0,to=6211050,by=47); samples[1]=1
samples=1:5000

p1=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                metric="sensitivity",y="Sensitivity") 
p2=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                metric="specificity",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                metric="precision",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,dist_list,
                metric="accuracy",y="accuracy")

##Saving "sensitivity" only
dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"sensitivity_diffFitness_first5000.pdf",sep="/"),plot=p1,device="pdf",width=3*4,height=2*4)



require(gridExtra)
p_all=grid.arrange(p1,p2,p3,p4, ncol=2) #When viewed in small window the dashed lines don't show properly

#dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
#ggsave(file=paste(dir_of_workingScript,"metrics_diffPhenotypeData_MI.pdf",sep="/"),plot=p_all,device="pdf",width=3*6,height=2*6)

##When samples=1:300 or =1:5000 it cuts the edges so I changed the width and hieight a little bit
#ggsave(file=paste(dir_of_workingScript,"metrics_diffPhenotypeData_MI_first3000.pdf",sep="/"),plot=p_all,device="pdf",width=3*7,height=2*7)
#ggsave(file=paste(dir_of_workingScript,"metrics_diffPhenotypeData_MI_first300.pdf",sep="/"),plot=p_all,device="pdf",width=3*7,height=2*7)



#plot the result in a different way:

#Here I use different no. of genes for different lines. Therefore I re-prepare the objs needed
#================================================================================================================================================
#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc","mi_ternary","mi_ternary_collapsedCond")]
TF=( rowSums(df[,-c(6,7,8)])==5 )
sum(TF) #106



Original_pcc=cumsum( TF[order(df$pcc)])
Original_pcc_confusionMatrix_metrics=confusionMatrix_metrics(Original_pcc)

#(!)Have to verify whether the mi, mi_ternary_collapsedCond used Ternary data with 324 cutoffs
mi_ternary=cumsum( TF[order(df$mi_ternary)]) 
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt with
mi_ternary_confusionMatrix_metrics=confusionMatrix_metrics(mi_ternary)


mi_ternary_collapsedCond=cumsum( TF[order(df$mi_ternary_collapsedCond)])
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt with
mi_ternary_collapsedCond_confusionMatrix_metrics=confusionMatrix_metrics(mi_ternary_collapsedCond)


#load("Data/NichDeut_mi_ternary_allAnnot.RData") #I used |fitness|>=3 for this experiment (refer to: deutchobaur_et_al/correlationVSannotation_test.R).
#NichDeut_mi_ternary_confusionMatrix_metrics=confusionMatrix_metrics(NichDeut_mi_ternary_allAnnot$cumsum)


dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
               "Ternary data - MI"=mi_ternary_confusionMatrix_metrics,
               "Ternary data with condition collapsed - MI"=mi_ternary_collapsedCond_confusionMatrix_metrics)

#dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
#                     "Ternary data - MI"=mi_ternary_confusionMatrix_metrics,
#                     "2 combined ternary data - MI"=NichDeut_mi_ternary_confusionMatrix_metrics) #note this guy has different no. of total pairs


#================================================================================================================================================



#Complete graph sampling from all 7914231 pairs
#============================================================================================================================================
metricByCutoff=function(cutoffIndices,metric_list){
  
  result_list=list()
  
  for(i in seq(metric_list)){
    
    result_list[[i]]=data.frame(sensitivity=metric_list[[i]]$sensitivity[cutoffIndices],
                                specificity=metric_list[[i]]$specificity[cutoffIndices],
                                precision=metric_list[[i]]$precision[cutoffIndices],
                                accuracy=metric_list[[i]]$accuracy[cutoffIndices],
                                group=names(metric_list)[i])
  }
  
  result=Reduce(rbind,result_list) 
  result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend
  
  return(result)
}



samples=seq(from=0,to=7914231,by=1989)
samples=samples[-1] #throw out the first index because it would cause problem when doing log() on plotting
samples=c(2,seq(from=5,to=1989,by=5),samples) #I don't want to miss points using stringent cutoffs (ranking of of pairs <1989)


cutoffIndices=samples
result=metricByCutoff(cutoffIndices=cutoffIndices,metric_list=dist_list)

Cols=c("#09D38A","#EC8D3A","#832B05")



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
  #scale_y_continuous(trans='log10',limits = c(0.001,1))+
  #scale_x_continuous(trans='log10',limits = c(0.001,1))+
  scale_color_manual(values=Cols)+
  labs(color="Type of data - Correlation")+
  coord_flip() #this flips x and y axes

#============================================================================================================================================


#A graph that only shows good cutoffs
#============================================================================================================================================
metricByCutoff_2=function(cutoffIndices,metric_list){
  
  result_list=list()
  
  for(i in seq(metric_list)){
    
    result_list[[i]]=data.frame(sensitivity=metric_list[[i]]$sensitivity[cutoffIndices[i]],
                                specificity=metric_list[[i]]$specificity[cutoffIndices[i]],
                                precision=metric_list[[i]]$precision[cutoffIndices[i]],
                                accuracy=metric_list[[i]]$accuracy[cutoffIndices[i]],
                                group=names(metric_list)[i])
  }
  
  result=Reduce(rbind,result_list) 
  result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend
  
  return(result)
}

optimal_indices=c()
for(i in seq(dist_list)){
  optimal_indices[i]=which(dist_list[[i]]$precision==max(dist_list[[i]]$precision))
}



#Find the point for Ternary MI where it almost overlap with the optimum of of Ternary MI - collapsed condition
precision=dist_list$`Ternary data - MI`$precision
sensitivity=dist_list$`Ternary data - MI`$sensitivity
abs_diffSensitivity=abs(sensitivity - dist_list$`Ternary data with condition collapsed - MI`$sensitivity[optimal_indices[3]])
Index=which(diffSensitivity==min(abs_diffSensitivity)) #3811,3812....3841. I think I can simply pick anyone. Eg. 3188
optimal_indices[2]=3811



result_2=metricByCutoff_2(cutoffIndices=optimal_indices,metric_list=dist_list)



Cols=c("#09D38A","#EC8D3A","#832B05")

p=ggplot(data=result_2,aes(sensitivity,precision,color=group))+
  theme_minimal()+
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

#============================================================================================================================================


currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
pdf(file=paste(currentDir,"sensitivity_precision_diffFitness.pdf",sep="/"),width=3*4,height=2*4)
p
dev.off()



#ROC curve (x: FPR = 1-specificity, y: TPR= sensitivity) (ref: https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
#What are good AUCs?: http://gim.unmc.edu/dxtests/roc3.htm


##Set no. of samples to be plotted
#samples=1:7914231
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1
#Price_samples=seq(from=0,to=6211050,by=47); samples[1]=1
samples=1:6211050

plot_ROC=function(samples,size,dist_list,legend_title="Type of data - Correlation"){
  
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#09D38A","#5D46E0","#832B05","#EC8D3A","red") 
  
  
  
  
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
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagnal line
    #scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    #scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)+
    #scale_x_continuous(expand = c(0,0),labels=comma)+
    scale_color_manual(values=Cols) 
}


roc_diffFitness=plot_ROC(samples,size=1,dist_list)

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)

start.time=Sys.time()
ggsave(file=paste(dir_of_workingScript,"roc_diffFitness_with_Price.png",sep="/"),plot=roc_diffFitness,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension
end.time=Sys.time()
end.time-start.time #Time difference of 3.322486 mins



#Calculate AUC (Area under curve) using a package
library(pROC) #ref: https://stackoverflow.com/questions/4903092/calculate-auc-in-r

df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")]
response=( rowSums(df)==5 )


#pcc
predictor=1-strain1strain2_allAnnotations_allDistances$pcc

start.time = Sys.time()
pccAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time #Time difference of 22.78993 secs

pccAUC #Area under the curve: 0.9131

#ternary MI
predictor=1-strain1strain2_allAnnotations_allDistances$mi_ternary

start.time = Sys.time()
ter_miAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 

ter_miAUC #Area under the curve: 0.8432

#ternary-collapsed conditions MI
predictor=1-strain1strain2_allAnnotations_allDistances$mi_ternary_collapsedCond

start.time = Sys.time()
ter_collapsed_miAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 

ter_collapsed_miAUC #Area under the curve: 0.8404

#ternary-collapsed conditions mhd3
predictor=1-strain1strain2_allAnnotations_allDistances$mhd3_collapsedCond

start.time = Sys.time()
mhd3_collapsedCondAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 


mhd3_collapsedCondAUC #Area under the curve: 0.8537





#ternary MI (2 combined datasets)
predictor=NichDeut_mi_ternary_allAnnot$mi_ternary
response=NichDeut_mi_ternary_allAnnot$sameORnot

start.time = Sys.time()
two_ter_miAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 

two_ter_miAUC #Area under the curve: 0.8757





#PR curve (x: recall(sensitivity), y: precision  ) (ref: https://classeval.wordpress.com/introduction/introduction-to-the-precision-recall-plot/)

#samples=1:7914231
#samples=seq(from=0,to=7914231,by=13); samples[1]=1 #7914231=3^2*13*17*23*173
#Price_samples=seq(from=0,to=6211050,by=47); samples[1]=1
Price_samples=1:6211050

plot_PRC=function(samples,size,dist_list,legend_title="Type of data - Correlation"){
  
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#09D38A","#5D46E0","#832B05","#EC8D3A","red") 
  
  
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(recall=dist_list[[1]]$sensitivity[samples],precision=dist_list[[1]]$precision[samples]),aes(recall,precision,color=names(dist_list[1])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=dist_list[[2]]$sensitivity[samples],precision=dist_list[[2]]$precision[samples]),aes(recall,precision,color=names(dist_list[2])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=dist_list[[3]]$sensitivity[samples],precision=dist_list[[3]]$precision[samples]),aes(recall,precision,color=names(dist_list[3])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=dist_list[[4]]$sensitivity[samples],precision=dist_list[[4]]$precision[samples]),aes(recall,precision,color=names(dist_list[4])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=dist_list[[5]]$sensitivity[samples],precision=dist_list[[5]]$precision[samples]),aes(recall,precision,color=names(dist_list[5])),size=size,alpha=0.7)+
    
    
    #this is the negative control
    geom_abline(slope=0,intercept=(dist_list[[1]]$TP[1] + dist_list[[1]]$FN[1])/dim(dist_list[[1]])[1],linetype="dashed")+ #intercept= Positive/(Negative + Positive) = Positive/Total = (TP + FN)/total
    
    
    #scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    #scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="Recall (sensitivity)",y="Precision",aesthetic='custom text',color=legend_title)+
    #scale_x_continuous(expand = c(0,0),labels=comma)+
    scale_color_manual(values=Cols) 
}


prc_diffFitness=plot_PRC(Price_samples,size=1,dist_list)

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)

start.time=Sys.time()
ggsave(file=paste(dir_of_workingScript,"prc_diffFitness_withPrice.png",sep="/"),plot=prc_diffFitness,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension
end.time=Sys.time()
end.time-start.time #Time difference of 8.260464 mins


#Part of the PR curve:




#install.packages("PRROC")
library(PRROC)

df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")]
response=( rowSums(df)==5 )


#pcc
predictor=1-strain1strain2_allAnnotations_allDistances$pcc

pr.curve(scores.class0 = predictor,weights.class0=as.numeric(response))
##Area under curve (Integral): 0.01240156 
##Area under curve (Davis & Goadrich):0.01240133 



#ternary MI
predictor=1-strain1strain2_allAnnotations_allDistances$mi_ternary

pr.curve(scores.class0 = predictor,weights.class0=as.numeric(response))
##Area under curve (Integral): 0.004295017 
##Area under curve (Davis & Goadrich): 0.004318629 


#ternary-collapsed conditions MI
predictor=1-strain1strain2_allAnnotations_allDistances$mi_ternary_collapsedCond

pr.curve(scores.class0 = predictor,weights.class0=as.numeric(response))
##Area under curve (Integral): 0.003340061 
##Area under curve (Davis & Goadrich): 0.003332748 


#random (Positive/Total) * 1 :
106/7914231 #1.34e-05




