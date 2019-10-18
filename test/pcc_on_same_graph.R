pheno_dat=strain1strain2_allAnnotations_allDistances

CoAnnotation_list=list()
CoAnnotation_list[[1]]=ifelse(pheno_dat[,"Pwy"]==1,T,F)
CoAnnotation_list[[2]]=ifelse(pheno_dat[,"pcomplex"]==1,T,F)
CoAnnotation_list[[3]]=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex")])==2,T,F)
CoAnnotation_list[[4]]=ifelse(pheno_dat[,"operon"]==1,T,F)
CoAnnotation_list[[5]]=ifelse(pheno_dat[,"regulator"]==1,T,F)
CoAnnotation_list[[6]]=ifelse(pheno_dat[,"kegg_modules"]==1,T,F)
CoAnnotation_list[[7]]=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5,T,F)


distance="pcc"

cumsums=lapply(CoAnnotation_list,FUN=function(CoAnnotation){
  cumsum(CoAnnotation[order(pheno_dat[,distance])])
})


str(cumsums)


slopeForNegativeControl=1/7914231
#slopeForPositiveControl= #Unfortunately positive controal for different annotation sets are different becuase the total no. of co-annotations are different
library(scales) #This is for label=comma (not using scientific notation and put commas) to work

##The corresponding plot for "pcc"

samples=1:6000
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1
alpha=0.5
size=1

# Ranked pairs VS. cummulative sum/total co-annotation
ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,cumSum_ratio=cumsums[[1]][samples]/sum(CoAnnotation_list[[1]])),aes(no,cumSum_ratio,color="EcoCyc Pathway"),size=size,alpha=alpha)+
  geom_line(data=data.frame(no=samples,cumSum_ratio=cumsums[[2]][samples]/sum(CoAnnotation_list[[2]])),aes(no,cumSum_ratio,color="EcoCyc Protein Complex"),size=size,alpha=alpha)+
  geom_line(data=data.frame(no=samples,cumSum_ratio=cumsums[[3]][samples]/sum(CoAnnotation_list[[3]])),aes(no,cumSum_ratio,color="Intersection of Pathway and Protein Complex"),size=size,alpha=alpha)+
  geom_line(data=data.frame(no=samples,cumSum_ratio=cumsums[[7]][samples]/sum(CoAnnotation_list[[7]])),aes(no,cumSum_ratio,color="Intersection of all 5 annotation sets"),size=size,alpha=alpha)+
  
  ##I tried to show legend for geom_abline but even the following failed
  ##https://groups.google.com/forum/#!topic/ggplot2/HtMK2ynltyw
  geom_abline(intercept=0,slope=slopeForNegativeControl,color="#D76E6E",linetype='dashed',size=1)+ #negative control
  labs(color="Annotation set")+
  labs(x = "Ranking of Distance",y="Cumulative No. of co-annotated pairs / Total no. of co-annotated pairs",aesthetic='custom text')+
  scale_x_continuous(expand = c(0,0),labels=comma)
  
  
  
  







