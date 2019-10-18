#1. sensitivity VS. precision  2. ROC-AUC using Wang's GO distances


df=strain1strain2_allAnnotations_allDistances

summary(df$Wang_BP)
summary(df$Wang_MF)
summary(df$Wang_CC)


#pick an arbitrary cutoff for "co-annotation"
sum(df$Wang_BP>=.8,na.rm=T)/dim(df)[1] #0.006305224
sum(df$Wang_MF>=.8,na.rm=T)/dim(df)[1] #0.01988886
sum(df$Wang_CC>=.8,na.rm=T)/dim(df)[1] #0.08559088



get_plots=function(cutoff=.8,semanticSimilarity,samples,pcc=strain1strain2_allAnnotations_allDistances$pcc){
  
  coannotated=ifelse(semanticSimilarity>=cutoff & !is.na(semanticSimilarity),1,0)
  
  confusion_metrics=confusionMatrix_metrics(coannotated[order(pcc)] %>% cumsum)
  
  result=data.frame(sensitivity=confusion_metrics$sensitivity[samples],
                    specificity=confusion_metrics$specificity[samples],
                    precision=confusion_metrics$precision[samples],
                    accuracy=confusion_metrics$accuracy[samples])
  
  ## sensitivity VS. precision
  p1_1=ggplot(data=result,aes(sensitivity,precision))+
    theme_minimal()+
    geom_line()+
    theme(plot.title= element_text(size = 20),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          axis.text.x = element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title=element_text(size=20))+
    labs(color="Type of data - Correlation")+
    coord_flip() #this flips x and y axes
  
  ## sensitivity VS. precision - loglog transformation
  p1_2=ggplot(data=result,aes(sensitivity,precision))+
    theme_minimal()+
    geom_line()+
    theme(plot.title= element_text(size = 20),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          axis.text.x = element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title=element_text(size=20))+
    scale_y_continuous(trans='log10',limits = c(0.000001,1))+
    scale_x_continuous(trans='log10',limits = c(0.000001,1))+
    labs(color="Type of data - Correlation")+
    coord_flip() #this flips x and y axes
  
  
  ## ROC-AUC
  p2=ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(FPR=1-result$specificity,sensitivity=result$sensitivity),aes(FPR,sensitivity))+
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagnal line
    scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text')
  
  
  return(list(p1_1,p1_2,p2))
}





samples=seq(from=0,to=7914231,by=1989)
samples=samples[-1] #throw out the first index because it would cause problem when doing log() on plotting
samples=c(2,seq(from=5,to=1989,by=5),samples) #I don't want to miss points using stringent cutoffs (ranking of of pairs <1989)


BP_plots=get_plots(cutoff=.8,df$Wang_BP,samples=samples)
MF_plots=get_plots(cutoff=.8,df$Wang_MF,samples=samples)
CC_plots=get_plots(cutoff=.8,df$Wang_CC,samples=samples)


BP_plots_2=get_plots(cutoff=.9,df$Wang_BP,samples=samples)
MF_plots_2=get_plots(cutoff=.9,df$Wang_MF,samples=samples)
CC_plots_2=get_plots(cutoff=.9,df$Wang_CC,samples=samples)



#Conclusion: Only BP_plots show significance

##optimal cutoff for sensitivity VS. precision
cutoff=.9
coannotated=ifelse(semanticSimilarity>=cutoff & !is.na(semanticSimilarity),1,0)
confusion_metrics=confusionMatrix_metrics(coannotated[order(pcc)] %>% cumsum)

Index=1946 #optimal cutoff determined by using all 5 annots (|PCC|=0.67)
confusion_metrics$precision[1946] #0.06834532
confusion_metrics$sensitivity[1946] #0.004087277



##Calculate AUC for BP's ROC

###cutoff=.8
cutoff=.8
response=ifelse(df$Wang_BP>=cutoff & !is.na(df$Wang_BP),1,0)
predictor=strain1strain2_allAnnotations_allDistances$pcc

start.time = Sys.time()
BP_AUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time #Time difference of 1.137108 mins

BP_AUC #Area under the curve: 0.5195


###cutoff=.9
cutoff=.9
response=ifelse(df$Wang_BP>=cutoff & !is.na(df$Wang_BP),1,0)
predictor=strain1strain2_allAnnotations_allDistances$pcc

start.time = Sys.time()
BP_AUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time 

BP_AUC #Area under the curve: 0.5204



#iterate over a list of cutoffs
cutoffs=seq(from=0.05,to=1,by=0.05)
predictor=strain1strain2_allAnnotations_allDistances$pcc

start.time = Sys.time()
for(cutoff in cutoffs){
  response=ifelse(df$Wang_BP>=cutoff & !is.na(df$Wang_BP),1,0)
  
  auc(response,predictor) %>% print
}
end.time = Sys.time()
end.time - start.time #Time difference of 24.30096 mins


##Result:
##Area under the curve: 0.503
##Area under the curve: 0.5041
##Area under the curve: 0.5047
##Area under the curve: 0.5051
##Area under the curve: 0.5065
##Area under the curve: 0.5083
##Area under the curve: 0.5093
##Area under the curve: 0.5107
##Area under the curve: 0.5121
##Area under the curve: 0.5145
##Area under the curve: 0.5164
##Area under the curve: 0.5187
##Area under the curve: 0.52
##Area under the curve: 0.5194
##Area under the curve: 0.519
##Area under the curve: 0.5195
##Area under the curve: 0.5195
##Area under the curve: 0.5204
##Area under the curve: 0.5191
##Area under the curve: 0.5192



