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
plot_metrics=function(samples,alpha,size,four_annot_list,metric,x_lab="high |PCC| -- ranked pairs -- low |PCC|",y,legend_title="Annotation set(s)"){
  
  ##Cols=c("#56B4E9","#F3518A","#BE70EA","#09D38A") 
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#09D38A","#56B4E9","#BE70EA","#F3518A")  
  
  random_metric=paste("random",metric,sep="_")
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[1]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[1])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[2]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[2])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[3]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[3])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[4]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[4])),size=size)+
    
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[1]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[1])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[2]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[2])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[3]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[3])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[4]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[4])),size=size,alpha=alpha,linetype = "dashed")+
    #geom_abline(slope=1/max(samples),intercept=0,linetype="dashed")+ # (for sensitivity use this) This is the diagnal line
    
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)+
    scale_color_manual(values=rep(Cols,2)) 
}
#========================================================================================================================


#pwy
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcc")]
Pwy=cumsum(df$Pwy[order(df$pcc)])
Pwy_confusionMatrix_metrics=confusionMatrix_metrics(Pwy)



#pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("pcomplex","pcc")]
pcomplex=cumsum(df$pcomplex[order(df$pcc)])
pcomplex_confusionMatrix_metrics=confusionMatrix_metrics(pcomplex)



#pwy+pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","pcc")]
TF=( (df$Pwy+df$pcomplex)==2 )
sum(TF) #230
PwyANDpcomplex=cumsum( TF[order(df$pcc)])
PwyANDpcomplex_confusionMatrix_metrics=confusionMatrix_metrics(PwyANDpcomplex)



#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
TF=( rowSums(df[,-6])==5 )
sum(TF) #106
All_annotSet=cumsum( TF[order(df$pcc)])
All_annotSet_confusionMatrix_metrics=confusionMatrix_metrics(All_annotSet)




four_annot_list=list("Same Pathways"=Pwy_confusionMatrix_metrics,
                           "Same protein complexes"=pcomplex_confusionMatrix_metrics,
                           "Same pathways & protein complexes"=PwyANDpcomplex_confusionMatrix_metrics,
                           "Same in all 5 sets"=All_annotSet_confusionMatrix_metrics)


##Set no. of samples to be plotted
##samples=1:200
##samples=1:4000
##samples=1:5000
samples=seq(from=0,to=7914231,by=1989); samples[1]=1

p1=plot_metrics(samples=samples,alpha=0.5,size=1,four_annot_list,metric="sensitivity",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="Sensitivity")
p2=plot_metrics(samples=samples,alpha=0.5,size=1,four_annot_list,metric="specificity",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,four_annot_list,metric="precision",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,four_annot_list,metric="accuracy",x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="accuracy")


##Saving "sensitivity" only
#dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
#ggsave(file=paste(dir_of_workingScript,"sensitivity_diffAnnotation_first4000.pdf",sep="/"),plot=p1,device="pdf",width=3*4,height=2*4)
#ggsave(file=paste(dir_of_workingScript,"sensitivity_diffAnnotation_first5000.pdf",sep="/"),plot=p1,device="pdf",width=3*4,height=2*4)
#ggsave(file=paste(dir_of_workingScript,"sensitivity_diffAnnotation.pdf",sep="/"),plot=p1,device="pdf",width=3*4,height=2*4)



require(gridExtra)
p_all=grid.arrange(p1, p2,p3,p4, ncol=2) #When viewed in small window the dashed lines don't show properly

#dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
#ggsave(file=paste(dir_of_workingScript,"metrics_diffAnnotation_first4000.pdf",sep="/"),plot=p_all,device="pdf",width=3*6,height=2*6)


# (!) The following doesn't allow p_all to be properly generated
#pdf(file="metrics_diffAnnot_first4000",width=12*2,height=8*2) 
#p_all
#dev.off()

#pdf(file="metrics_diffAnnot",width=12,height=8)
#p_all
#dev.off()



#plot the result in a different way:

#sensitivity VS. precision (axes can be changed to have different metrics on x and y)

##Using various cutoffs
#============================================================================================================================================
metricByCutoff=function(cutoffs=1-c(0.9,0.8,0.7,0.6,0.1),metric_list,distance=pcc_dist){
  
  
  
  result_list=list()
  count=1
  for(i in seq(four_annot_list)){
    for(cutoff in cutoffs){
      
      Index=sum(pcc_dist<cutoff|pcc_dist==cutoff)
      result_list[[count]]=data.frame(sensitivity=four_annot_list[[i]]$sensitivity[Index],
                                      specificity=four_annot_list[[i]]$specificity[Index],
                                      precision=four_annot_list[[i]]$precision[Index],
                                      accuracy=four_annot_list[[i]]$accuracy[Index],
                                      group=names(four_annot_list)[i])
      
      count=count+1
    }
    
  }
  
  result=Reduce(rbind,result_list) 
  result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend
  
  return(result)
}

pcc_dist=sort(strain1strain2_allAnnotations_allDistances$pcc)


cutoffs=1-c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
result=metricByCutoff(metric_list=four_annot_list,cutoffs=cutoffs)

Cols=c("#56B4E9","#F3518A","#BE70EA","#09D38A") 


p=ggplot(data=result,aes(sensitivity,precision,color=group))+ #simply swaping sensitivity & precision will give weird plot (segments for the plot will be messed up because they link from samllest x value to largest x value)
  theme_minimal()+
  geom_line()+
  geom_point(size=3)+
  geom_text(label=c("0.9","","","0.6","","","","","0.1") %>% rep(4),
            vjust=rep(-1,4*length(cutoffs)),size=5,show.legend = FALSE)+ #show.legend = FALSE removes the mysterious "a" in the legend (https://stackoverflow.com/questions/18337653/remove-a-from-legend-when-using-aesthetics-and-geom-text)
  theme(plot.title= element_text(size = 20),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=20))+
  scale_y_continuous(trans='log10',limits = c(0.00001,1))+
  scale_x_continuous(trans='log10',limits = c(0.00001,1))+
  scale_color_manual(values=Cols)+
  labs(color="Annotation set(s)")+
  coord_flip() #this flips x and y axes


currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
pdf(file=paste(currentDir,"sensitivity_precision_Annot_multipleCutoff.pdf",sep="/"),width=3*4,height=2*4)
p
dev.off()
#============================================================================================================================================


##Using an optimum cutoff determined by using all 5 annot (cutoff: |PCC|=0.67 determined in confusioinMatrix_analysis_diffDistances.R)
#============================================================================================================================================

metricByCutoff=function(cutoffIndices,metric_list,distance=pcc_dist){
  
  result_list=list()
  
  for(i in seq(four_annot_list)){
    
    result_list[[i]]=data.frame(sensitivity=four_annot_list[[i]]$sensitivity[cutoffIndices],
                                specificity=four_annot_list[[i]]$specificity[cutoffIndices],
                                precision=four_annot_list[[i]]$precision[cutoffIndices],
                                accuracy=four_annot_list[[i]]$accuracy[cutoffIndices],
                                group=names(four_annot_list)[i])
  }
  
  result=Reduce(rbind,result_list) 
  result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend
  
  return(result)
}


pcc_dist=sort(strain1strain2_allAnnotations_allDistances$pcc)
#samples=seq(from=0,to=7914231,by=1989)
#samples=samples[-1] #throw out the first index because it would cause problem when doing log() on plotting
#samples=c(2,seq(from=5,to=1989,by=5),samples) #I don't want to miss points using stringent cutoffs (ranking of of pairs <1989)
samples=1946 #optimal cutoff determined by using all 5 annots (|PCC|=0.67)


cutoffIndices=samples

result=metricByCutoff(cutoffIndices=cutoffIndices,metric_list=four_annot_list)


#sensitivity VS. precision (axes can be changed to have different metrics on x and y)

Cols=c("#56B4E9","#F3518A","#BE70EA","#09D38A") 


p=ggplot(data=result,aes(sensitivity,precision,color=group))+ #simply swaping sensitivity & precision will give weird plot (segments for the plot will be messed up because they link from samllest x value to largest x value)
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
  labs(color="Annotation set(s)")+
  coord_flip() #this flips x and y axes


currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
pdf(file=paste(currentDir,"sensitivity_precision_Annot.pdf",sep="/"),width=3*4,height=2*4)
p
dev.off()

#============================================================================================================================================


#ROC curve (x: FPR = 1-specificity, y: TPR= sensitivity) (ref: https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
#What are good AUCs?: http://gim.unmc.edu/dxtests/roc3.htm

plot_ROC=function(samples,alpha,size,four_annot_list,legend_title="Annotation set(s)"){
  
  ##Cols=c("#56B4E9","#F3518A","#BE70EA","#09D38A") 
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#09D38A","#56B4E9","#BE70EA","#F3518A")  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(FPR=1-four_annot_list[[1]]$specificity[samples],sensitivity=four_annot_list[[1]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(four_annot_list[1])),size=size)+
    geom_line(data=data.frame(FPR=1-four_annot_list[[2]]$specificity[samples],sensitivity=four_annot_list[[2]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(four_annot_list[2])),size=size)+
    geom_line(data=data.frame(FPR=1-four_annot_list[[3]]$specificity[samples],sensitivity=four_annot_list[[3]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(four_annot_list[3])),size=size)+
    geom_line(data=data.frame(FPR=1-four_annot_list[[4]]$specificity[samples],sensitivity=four_annot_list[[4]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(four_annot_list[4])),size=size)+
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagnal line
    scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)+
    scale_color_manual(values=Cols) 
}


samples=seq(from=0,to=7914231,by=1989); samples[1]=1
roc_diffAnnot=plot_ROC(samples=samples,alpha=0.5,size=1,four_annot_list)

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"roc_diffAnnot.png",sep="/"),plot=roc_diffAnnot,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension

#Calculate AUC (Area under curve) using a package
library(pROC) #ref: https://stackoverflow.com/questions/4903092/calculate-auc-in-r

predictor=strain1strain2_allAnnotations_allDistances$pcc

#pwy
response=strain1strain2_allAnnotations_allDistances$Pwy

start.time = Sys.time()
pwyAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time #Time difference of 29.62732 secs

pwyAUC #Area under the curve: 0.5954


#pcomplex
response=strain1strain2_allAnnotations_allDistances$pcomplex

start.time = Sys.time()
pcomplexAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time #Time difference of 29.51495 secs

pcomplexAUC #Area under the curve: 0.6832


#pwy+pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","pcc")]
response=( (df$Pwy+df$pcomplex)==2 )

start.time = Sys.time()
pwyANDpcomAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time #Time difference of 22.95037 secs

pwyANDpcomAUC #Area under the curve: 0.7956



#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
response=( rowSums(df[,-6])==5 )

start.time = Sys.time()
allAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time #Time difference of 22.89565 secs

allAUC #Area under the curve: 0.9131




#Determination of the optimal cutoff
dat=four_annot_list$`Same in all 5 sets`
str(dat)

ROC_x=1-dat$specificity #1-specificity
ROC_y=dat$sensitivity #sensitivity

samples=seq(from=0,to=7914231,by=1989); samples[1]=1
plot(ROC_x[samples],ROC_y[samples])
abline(v=ROC_x[which.max((order(ROC_x,decreasing = F)+order(ROC_y,decreasing = F)))/2]) #optimal cutoff


optimalIndex=which.max((order(ROC_x,decreasing = F)+order(ROC_y,decreasing = F)))/2

sort(1-strain1strain2_allAnnotations_allDistances$pcc,decreasing = T)[optimalIndex] #at the optimal cutoff |pcc|=0.076 (seems odd)




