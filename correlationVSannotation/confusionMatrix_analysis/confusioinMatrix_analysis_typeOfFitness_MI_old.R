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







#Function to get the graph. 
#========================================================================================================================
plot_metrics=function(samples,alpha,size,three_dist_list,metric,x_lab="Low distance -- ranked pairs -- High distance",y,legend_title="Type of data - Correlation"){
  
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
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[1]][[metric]][samples]),aes(no,Metric,color=names(three_dist_list[1])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[2]][[metric]][samples]),aes(no,Metric,color=names(three_dist_list[2])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[3]][[metric]][samples]),aes(no,Metric,color=names(three_dist_list[3])),size=size)+
    
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[1]][[random_metric]][samples]),aes(no,Metric,color=names(three_dist_list[1])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[2]][[random_metric]][samples]),aes(no,Metric,color=names(three_dist_list[2])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[3]][[random_metric]][samples]),aes(no,Metric,color=names(three_dist_list[3])),size=size,alpha=alpha,linetype = "dashed")+
    
    
    
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


mi_ternary=cumsum( TF[order(df$mi_ternary)]) 
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt with
mi_ternary_confusionMatrix_metrics=confusionMatrix_metrics(mi_ternary,total=6211050)


mi_ternary_collapsedCond=cumsum( TF[order(df$mi_ternary_collapsedCond)])
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt with
mi_ternary_collapsedCond_confusionMatrix_metrics=confusionMatrix_metrics(mi_ternary_collapsedCond,total=6211050)


load("Data/NichDeut_mi_ternary_allAnnot.RData") #I used |fitness|>=3 for this experiment (refer to: deutchobaur_et_al/correlationVSannotation_test.R).
NichDeut_mi_ternary_confusionMatrix_metrics=confusionMatrix_metrics(NichDeut_mi_ternary_allAnnot$cumsum,total=6211050)


#three_dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
#                     "Ternary data - MI"=mi_ternary_confusionMatrix_metrics,
#                     "Ternary data with condition collapsed - MI"=mi_ternary_collapsedCond_confusionMatrix_metrics)

three_dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
                     "Ternary data - MI"=mi_ternary_confusionMatrix_metrics,
                     "2 combined ternary data - MI"=NichDeut_mi_ternary_confusionMatrix_metrics)

##Set no. of samples to be plotted
##samples=1:300
##samples=1:4000
##samples=seq(from=0,to=7914231,by=1989); samples[1]=1
#samples=seq(from=0,to=6211050,by=47); samples[1]=1
samples=1:5000

p1=plot_metrics(samples=samples,alpha=0.5,size=1,three_dist_list,
                metric="sensitivity",y="Sensitivity") 
p2=plot_metrics(samples=samples,alpha=0.5,size=1,three_dist_list,
                metric="specificity",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,three_dist_list,
                metric="precision",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,three_dist_list,
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
metricByCutoff=function(cutoffs_list,metric_list,distance_list){
  
  
  
  result_list=list()
  count=1
  for(i in seq(metric_list)){
    cutoffs=cutoffs_list[[i]]
    
    for(cutoff in cutoffs){
      Index=sum(distance_list[[i]]<cutoff|distance_list[[i]]==cutoff)
      
      result_list[[count]]=data.frame(sensitivity=metric_list[[i]]$sensitivity[Index],
                                      specificity=metric_list[[i]]$specificity[Index],
                                      precision=metric_list[[i]]$precision[Index],
                                      accuracy=metric_list[[i]]$accuracy[Index],
                                      group=names(metric_list)[i])
      
      count=count+1
    }
    
  }
  
  result=Reduce(rbind,result_list) 
  result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend
  
  return(result)
}


#parameter prep
#============================================================================================================================================
dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
               "Ternary data - MI"=mi_ternary_confusionMatrix_metrics,
               "2 combined ternary data - MI"=NichDeut_mi_ternary_confusionMatrix_metrics #note this guy has different no. of total pairs
               )

pcc_dist=sort(strain1strain2_allAnnotations_allDistances$pcc)
mi_ternary_dist=sort(strain1strain2_allAnnotations_allDistances$mi_ternary)
mi_2ternary_dist=sort(NichDeut_mi_ternary_allAnnot$mi_ternary)


distance_list=list(pcc_dist,mi_ternary_dist,mi_2ternary_dist)

cutoffs_list=list(
  1-c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1), #pcc #doing 1-c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1) explains why we choose pcc=0.6 as the cutoff
  1-c(0.2,0.1,0.01,0.001), #mi for tenary (# max(mi)=0.7236666, so I can't use 1-0.9 as the first cutoff )
  1-c(0.2,0.1,0.01,0.001) #mi for 2 ternary (# max(mi)=0.7333098, so I can't use 1-0.9 as the first cutoff )
)
#============================================================================================================================================


result=metricByCutoff(cutoffs_list=cutoffs_list,metric_list=dist_list,distance_list=distance_list)


#sensitivity VS. precision (axes can be changed to have different metrics on x and y)
n_cutoffs=length(unlist(cutoffs_list))
Cols=c("#09D38A","#EC8D3A","#832B05")


p=ggplot(data=result,aes(sensitivity,precision,color=group))+
  theme_minimal()+
  geom_line()+
  geom_point(size=3)+
  geom_text(label=c(c("0.9","","","0.6","","","","","0.1"),
                               c("0.2","0.1","","0.001"),
                               c("0.2","0.1","","0.001")),
            vjust=rep(-1,n_cutoffs),size=5,show.legend = FALSE)+ #show.legend = FALSE removes the mysterious "a" in the legend (https://stackoverflow.com/questions/18337653/remove-a-from-legend-when-using-aesthetics-and-geom-text)
  theme(plot.title= element_text(size = 20),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=20))+
  scale_y_continuous(trans='log10',limits = c(0.00001,1))+
  scale_x_continuous(trans='log10',limits = c(0.00001,1))+
  scale_color_manual(values=Cols)+
  labs(color="Type of data - Correlation")+
  coord_flip() #this flips x and y axes
  



#currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
pdf(file=paste(currentDir,"sensitivity_precision_combinedFitness.pdf",sep="/"),width=3*4,height=2*4)
p
dev.off()



#ROC curve (x: FPR = 1-specificity, y: TPR= sensitivity) (ref: https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
#What are good AUCs?: http://gim.unmc.edu/dxtests/roc3.htm


#Can I use different no. of genes for different lines?

plot_ROC=function(samples,alpha,size,three_dist_list,legend_title="Type of data - Correlation"){
  
  #Cols=c("#09D38A","#EC8D3A","#832B05")
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#832B05","#09D38A","#EC8D3A") 
  
  
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(FPR=1-three_dist_list[[1]]$specificity[samples],sensitivity=three_dist_list[[1]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(three_dist_list[1])),size=size)+
    geom_line(data=data.frame(FPR=1-three_dist_list[[2]]$specificity[samples],sensitivity=three_dist_list[[2]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(three_dist_list[2])),size=size)+
    geom_line(data=data.frame(FPR=1-three_dist_list[[3]]$specificity[samples],sensitivity=three_dist_list[[3]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(three_dist_list[3])),size=size)+
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagnal line
    
    
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)+
    scale_color_manual(values=Cols) 
}

samples=seq(from=0,to=6211050,by=47); samples[1]=1
roc_diffFitness=plot_ROC(samples=samples,alpha=0.5,size=1,three_dist_list)

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"roc_diffFitness.png",sep="/"),plot=roc_diffFitness,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension


#Calculate AUC (Area under curve) using a package
library(pROC) #ref: https://stackoverflow.com/questions/4903092/calculate-auc-in-r

df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
response=( rowSums(df[,-6])==5 )


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


#ternary-collapsed conditions MI
predictor=strain1strain2_allAnnotations_allDistances$mi_ternary_collapsedCond

start.time = Sys.time()
terCollapsed_miAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 

terCollapsed_miAUC #Area under the curve: 0.8404




#ternary mhd3
predictor=strain1strain2_allAnnotations_allDistances$mhd3

start.time = Sys.time()
mhd3AUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 


mhd3AUC #Area under the curve: 0.8546


#ternary-collapsed conditions mhd3
predictor=strain1strain2_allAnnotations_allDistances$mhd3_collapsedCond

start.time = Sys.time()
mhd3_collapsedCondAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 


mhd3_collapsedCondAUC #Area under the curve: 0.8537













