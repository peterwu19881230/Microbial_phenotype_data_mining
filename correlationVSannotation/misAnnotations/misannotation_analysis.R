#analyze misannotation sets

load("Data/strain1_strain2_misCoannotation.RData")
str(strain1_strain2_misCoannotation)


#Here I will only use intersection of all annotation sets. I might be doing others if needed.

##calculate cumulative sum of "in all 5 annotation sets"
load("Data/strain1strain2_allDistances.RData")
str(strain1strain2_allDistances)


strain1_strain2_misCoannotation_allDistance=merge(strain1_strain2_misCoannotation,strain1strain2_allDistances ,by=c("strain1","strain2"))
str(strain1_strain2_misCoannotation_allDistance)


##calculate metrics based on confusion matrix
percentMis=c("0%", "0.5%", "5%", "25%", "50%")
type=c("shuff","remov")

df=strain1_strain2_misCoannotation_allDistance

start.time=Sys.time()

col_names=c()
for(i in 1:5){
  for(j in 1:2){
    coAnnotColumns=df[,grepl(paste("_",percentMis[i],"_",sep=""),names(df)) & grepl(type[j],names(df))]
    coAnnotColumn=ifelse(rowSums(coAnnotColumns)==5,1,0)
    
    
    if(i==1 & j==1){
      inAll5_column=coAnnotColumn
    }else{
      inAll5_column=cbind(inAll5_column,coAnnotColumn)
    }
    
    col_names=c(col_names,paste(c(percentMis[i],type[j]),collapse="_"))
    
  }
  
}

colnames(inAll5_column)=col_names


end.time=Sys.time()
end.time-start.time #Time difference of 2.547428 mins


inAll5_column_ordered=inAll5_column[order(strain1_strain2_misCoannotation_allDistance$pcc),] 

apply(inAll5_column_ordered,2,sum)









##analysis by ROC 
cumsums=apply(inAll5_column_ordered,2,cumsum)
conf_metrics_list=apply(cumsums,2,confusionMatrix_metrics)


plot_ROC=function(samples,alpha,size,conf_metrics_list,legend_title="Mis-annotation set(s)"){
  
  
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
    
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagonal line
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)
    
    return(p)
}


samples=seq(from=0,to=7914231,by=1989); samples[1]=1
roc_shuffle=plot_ROC(samples=samples,alpha=0.8,size=1,conf_metrics_list[grepl("shuff",names(conf_metrics_list))])
roc_remove=plot_ROC(samples=samples,alpha=0.8,size=1,conf_metrics_list[grepl("remov",names(conf_metrics_list))])

roc_shuffle
roc_remove




#Calculate AUC (Area under curve) using a package
library(pROC) #ref: https://stackoverflow.com/questions/4903092/calculate-auc-in-r



apply(inAll5_column,2,sum)  


predictor=strain1_strain2_misCoannotation_allDistance$pcc
for(i in (1:dim(inAll5_column)[2])[ apply(inAll5_column,2,sum)!=0 ]){  #only apply to columns that have at least one 1
  
  response=inAll5_column[,i]
  
  start.time = Sys.time()
  
  colnames(inAll5_column)[i] %>% print
  auc(response,predictor) %>% print #Syntax: auc(response,predictor)
  
  end.time = Sys.time()
  end.time - start.time #Time difference of 22.89565 secs
  
}

apply(inAll5_column,2,sum)


#no. of coannotations (intersection of 5 annotation sets)
#0%_shuff   0%_remov 0.5%_shuff 0.5%_remov   5%_shuff   5%_remov  25%_shuff  25%_remov 50%_shuff  50%_remov
#106        106        106        105         69        101          9         18          2          3
 

#Area under the curve 

#0%_shuff  0%_remov 0.5%_shuff  0.5%_remov   5%_shuff   5%_remov  25%_shuff  25%_remov 50%_shuff  50%_remov
# 0.9131     0.9131    0.9131     0.9127     0.9291      0.9251    0.8982      0.9862     0.8217    0.9998



#compare 0%_shuff (original) with 25%_shuff

##look at the |pcc| of the co-annotated pairs
pcc=1-strain1_strain2_misCoannotation_allDistance$pcc

pcc_original=pcc[(inAll5_column[,"0%_shuff"] %>% as.logical)]
pcc_25percent=pcc[(inAll5_column[,"25%_shuff"] %>% as.logical)]

hist(pcc_original, col=rgb(1,0,0,0.5), main="Overlapping Histogram", xlab="pcc")
hist(pcc_25percent, col=rgb(0,0,1,0.5), add=T)








