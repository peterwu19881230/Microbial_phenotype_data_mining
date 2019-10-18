
#Load data and define parameters
#===============================
load("Data/Price_pairwisePCC_longTable_coannotations.RData") #the pairwise price data here contains only strains that are also in Nichols
load("Data/NichDeut_mi_ternary.RData") #ternary_mi here means the mi based distance

NichDeut_mi_ternary$strain1=as.numeric(NichDeut_mi_ternary$strain1)
NichDeut_mi_ternary$strain2=as.numeric(NichDeut_mi_ternary$strain2)
temp=dplyr::left_join(NichDeut_mi_ternary,Price_pairwisePCC_longTable_coannotations[,-3],by=c("strain1","strain2"))
##apply(temp,2,anyNA) #there are NAs. I have to swap strain1 and strain2 for proper joining
temp[is.na(temp$Pwy),c("strain1","strain2")]=temp[is.na(temp$Pwy),c("strain2","strain1")] #swap strain1, strain2
temp=left_join(temp[,c("strain1","strain2","ternary_mi")],Price_pairwisePCC_longTable_coannotations[,-3],by=c("strain1","strain2"))
##apply(temp,2,anyNA) #all NAs are gone

large_df=temp[order(temp$ternary_mi),]
distance="ternary_mi"
filename_partial="Nich_Price"
samples=seq(from=0,to=6211050,by=47); samples[1]=1
text_size=25
#===============================


#(using pcc)function to get metrics from confusion matrix
get_metrics=function(complete_table,annot){
  df=complete_table[,c(annot,distance)]
  cumsum_=cumsum(df[,annot][order(df[,distance])])
  result=confusionMatrix_metrics(cumsum_)
  
  return(result)
}


annot_list=list()
for(annot in c("Pwy","pcomplex","operon","regulator","kegg_modules")){
  
  annot_list[[annot]]=get_metrics(large_df,annot)
}


#append "any coannotation" cumsum
df=large_df[,c("Pwy","pcomplex","operon","regulator","kegg_modules",distance)]
TF=( rowSums(df[,-6])>=1 )
##sum(TF) #223165
any_annotSet=cumsum( TF[order(df[,distance])])
any_annotSet_confusionMatrix_metrics=confusionMatrix_metrics(any_annotSet)

annot_list=c(annot_list,any=list(any_annotSet_confusionMatrix_metrics))



#Function to get the graph. 
#========================================================================================================================
plot_metrics=function(samples,alpha,size,four_annot_list,metric,x_lab="high similarity -- ranked pairs -- low similarity",y,legend_title="Annotation set(s)"){
  
  random_metric=paste("random",metric,sep="_")
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=text_size),
                                 legend.text=element_text(size=text_size),
                                 axis.text.x = element_text(size=text_size),
                                 axis.text.y=element_text(size=text_size),
                                 axis.title=element_text(size=text_size))+ 
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[1]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[1])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[2]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[2])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[3]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[3])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[4]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[4])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[5]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[5])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[6]][[metric]][samples]),aes(no,Metric,color=names(four_annot_list[6])),size=size)+
    
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[1]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[1])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[2]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[2])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[3]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[3])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[4]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[4])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[5]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[5])),size=size,alpha=alpha,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[6]][[random_metric]][samples]),aes(no,Metric,color=names(four_annot_list[6])),size=size,alpha=alpha,linetype = "dashed")+
    
    
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)
}
#========================================================================================================================





p1=plot_metrics(samples=samples,alpha=0.5,size=1,annot_list,metric="sensitivity",x_lab="high similarity -- ranked pairs -- low similarity",y="Sensitivity")
p2=plot_metrics(samples=samples,alpha=0.5,size=1,annot_list,metric="precision",x_lab="high similarity -- ranked pairs -- low similarity",y="precision")

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"/",filename_partial,"_ranking_sensitivity.png",sep=""),plot=p1,device="png",width=3*4,height=2*4)
ggsave(file=paste(dir_of_workingScript,"/",filename_partial,"_ranking_precision.png",sep=""),plot=p2,device="png",width=3*4,height=2*4)



p1_subset=plot_metrics(samples=1:5000,alpha=0.5,size=1,annot_list,metric="sensitivity",x_lab="high similarity -- ranked pairs -- low similarity",y="Sensitivity")
p2_subset=plot_metrics(samples=1:500,alpha=0.5,size=1,annot_list,metric="precision",x_lab="high similarity -- ranked pairs -- low similarity",y="precision")

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"/",filename_partial,"_ranking_sensitivity_subset.png",sep=""),plot=p1_subset,device="png",width=3*4,height=2*4)
ggsave(file=paste(dir_of_workingScript,"/",filename_partial,"_ranking_precision_subset.png",sep=""),plot=p2_subset,device="png",width=3*4,height=2*4)


#PR curve
plot_PRC=function(samples,size,annot_list,legend_title="Type of data - Correlation"){
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=text_size),
                                 legend.text=element_text(size=text_size),
                                 axis.text.x = element_text(size=text_size),
                                 axis.text.y=element_text(size=text_size),
                                 axis.title=element_text(size=text_size))+ 
    geom_line(data=data.frame(recall=annot_list[[1]]$sensitivity[samples],precision=annot_list[[1]]$precision[samples]),aes(recall,precision,color=names(annot_list[1])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=annot_list[[2]]$sensitivity[samples],precision=annot_list[[2]]$precision[samples]),aes(recall,precision,color=names(annot_list[2])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=annot_list[[3]]$sensitivity[samples],precision=annot_list[[3]]$precision[samples]),aes(recall,precision,color=names(annot_list[3])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=annot_list[[4]]$sensitivity[samples],precision=annot_list[[4]]$precision[samples]),aes(recall,precision,color=names(annot_list[4])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=annot_list[[5]]$sensitivity[samples],precision=annot_list[[5]]$precision[samples]),aes(recall,precision,color=names(annot_list[5])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=annot_list[[6]]$sensitivity[samples],precision=annot_list[[6]]$precision[samples]),aes(recall,precision,color=names(annot_list[6])),size=size,alpha=0.7)+
    
    
    #this is the negative control
    geom_abline(slope=0,intercept=(annot_list[[1]]$TP[1] + annot_list[[1]]$FN[1])/dim(annot_list[[1]])[1],linetype="dashed")+ #intercept= Positive/(Negative + Positive) = Positive/Total = (TP + FN)/total
    
    
    #scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    #scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="Recall (sensitivity)",y="Precision",aesthetic='custom text',color=legend_title)
  #+scale_x_continuous(expand = c(0,0),labels=comma)
  
}



prc_diffAnnot=plot_PRC(samples=samples,size=1,annot_list)
prc_diffAnnot

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"/",filename_partial,"_prc_diffAnnot.png",sep=""),plot=prc_diffAnnot,device="png",width=3*4,height=2*4)




#samples=1:5000
#prc_diffAnnot_subset=plot_PRC(samples,size=1,annot_list)
#prc_diffAnnot_subset






