#gene in pwy only

#====================================================================================
dat=strain1strain2_allAnnotations_allDistances
idInPwyCoannot=c(dat[dat$Pwy==1,"strain1"],dat[dat$Pwy==1,"strain2"]) %>% unique
length(idInPwyCoannot) #877

dat2=dat[ (dat$strain1 %in% idInPwyCoannot) & (dat$strain2 %in% idInPwyCoannot),]

coAnnotInPwy=cumsum(dat2$Pwy[order(dat2$pcc)])

coAnnotInPwy_confusionMatrix_metrics=confusionMatrix_metrics(coAnnotInPwy)

str(coAnnotInPwy_confusionMatrix_metrics) #384126 obs. of  16 variables
#====================================================================================

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


#pwy or pwcomplex or KEGG modules or regulon or operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
TF=( rowSums(df[,-6])>=1 )
sum(TF) #223165
ANY_annotSet=cumsum( TF[order(df$pcc)])
ANY_annotSet_confusionMatrix_metrics=confusionMatrix_metrics(ANY_annotSet)

annot_list=list("Same Pathways"=Pwy_confusionMatrix_metrics,
                "Same protein complexes"=pcomplex_confusionMatrix_metrics,
                "Same pathways & protein complexes"=PwyANDpcomplex_confusionMatrix_metrics,
                "Same in all 5 sets"=All_annotSet_confusionMatrix_metrics,
                "Same in any of 5 sets"=ANY_annotSet_confusionMatrix_metrics,
                "Same Pathways using pwy-coAnnotated pairs"=coAnnotInPwy_confusionMatrix_metrics)



plot_ROC=function(samples,alpha,size,annot_list,legend_title="Annotation set(s)"){
  
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#09D38A","grey","#56B4E9","#BE70EA","black","#F3518A")  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  
  p=ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+
    geom_line(data=data.frame(FPR=1-annot_list[[1]]$specificity[samples],sensitivity=annot_list[[1]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(annot_list[1])),size=size)+
    geom_line(data=data.frame(FPR=1-annot_list[[2]]$specificity[samples],sensitivity=annot_list[[2]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(annot_list[2])),size=size)+
    geom_line(data=data.frame(FPR=1-annot_list[[3]]$specificity[samples],sensitivity=annot_list[[3]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(annot_list[3])),size=size)+
    geom_line(data=data.frame(FPR=1-annot_list[[4]]$specificity[samples],sensitivity=annot_list[[4]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(annot_list[4])),size=size)+
    geom_line(data=data.frame(FPR=1-annot_list[[5]]$specificity[samples],sensitivity=annot_list[[5]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(annot_list[5])),size=size)+
    
    geom_line(data=data.frame(FPR=1-annot_list[[6]]$specificity,sensitivity=annot_list[[6]]$sensitivity),aes(FPR,sensitivity,color=names(annot_list[6])),size=size)+
    
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagnal line
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)+
    scale_color_manual(values=Cols) 
    
    return(p)
}

samples=seq(from=0,to=7914231,by=1989); samples[1]=1
plot_ROC(samples=samples,alpha=0.5,size=1,annot_list)



