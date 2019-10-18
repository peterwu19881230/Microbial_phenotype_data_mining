##change pheno_dat to get 2 different results (1 is using all the strains and the other is using strains that have at least 1 significant phenotype)
pheno_dat=strain1strain2_allAnnotations_allDistances
#pheno_dat=sigPheno_strain1strain2_allAnnotations_allDistances


#AnyCoAnnotation or other combinations of databases (later the experiment can be expanded to include any 2,3 databases)

CoAnnotation_dataframe=data.frame(any=rep(0,7914231),Pwy=rep(0,7914231),pcomplex=rep(0,7914231),operon=rep(0,7914231),regulator=rep(0,7914231),kegg_modules=rep(0,7914231),Pwy_and_kegg_modules=rep(0,7914231),
                                  Pwy_or_kegg_modules=rep(0,7914231),anyTwo=rep(0,7914231),anyThree=rep(0,7914231),anyFour=rep(0,7914231),anyFive=rep(0,7914231))
CoAnnotation_dataframe$any=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>0,T,F)
CoAnnotation_dataframe$Pwy=ifelse(pheno_dat[,"Pwy"]>0,T,F)
CoAnnotation_dataframe$pcomplex=ifelse(pheno_dat[,"pcomplex"]>0,T,F)
CoAnnotation_dataframe$operon=ifelse(pheno_dat[,"operon"]>0,T,F)
CoAnnotation_dataframe$regulator=ifelse(pheno_dat[,"regulator"]>0,T,F)
CoAnnotation_dataframe$kegg_modules=ifelse(pheno_dat[,"kegg_modules"]>0,T,F)
CoAnnotation_dataframe$Pwy_and_kegg_modules=ifelse(rowSums(pheno_dat[,c("Pwy","kegg_modules")])==2,T,F)
CoAnnotation_dataframe$Pwy_or_kegg_modules=ifelse(rowSums(pheno_dat[,c("Pwy","kegg_modules")])>0,T,F)

CoAnnotation_dataframe$anyTwo=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>=2,T,F)
CoAnnotation_dataframe$anyThree=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>=3,T,F)
CoAnnotation_dataframe$anyFour=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>=4,T,F)
CoAnnotation_dataframe$anyFive=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5,T,F)




setwd("/Users/peterwu/Dropbox/Nichols_Data_mining/correlationVSannotation/figs")

#Print every fig to a file
distances=c("pcc","spearman","euclidean","euclidean_qualitative","euclidean_collapsedCond","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond")

  
#reorder by ranking of PCC and calculate cumsum
enrich_CoAnnotation_list=list()

count=1
for(distance in distances){

  CoAnnotation=CoAnnotation_dataframe[order(pheno_dat[,distance]),]
  #Before I had: Error: vector memory exhausted (limit reached?)
  ##Solved by creating a new file: 
  ##https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
  
  
  ##Have to convert the logical TF to numeric TF in order for cumsum() to work
  CoAnnotation=convertMatOrDF(CoAnnotation,as.numeric)
  
  ##Have to divide cumsum() by randomExp for every annotation set (or combinations of sets)
  RandomExp=matrix(nrow=dim(CoAnnotation_dataframe)[1],ncol=dim(CoAnnotation_dataframe)[2])
  for(i in 1:dim(CoAnnotation_dataframe)[2]){
    RandomExp[,i]=sum(CoAnnotation_dataframe[,i])/length(CoAnnotation_dataframe[,i])*1:length(CoAnnotation_dataframe[,i])
  }
  RandomExp=as.data.frame(RandomExp)
  
  #pairwise divsion (not matrix division)
  enrich_CoAnnotation_list[[count]]=cumsum(CoAnnotation)/RandomExp
  
  count=count+1
}
names(enrich_CoAnnotation_list)=distances




#plot it
samples=1:4000
##samples=seq(from=0,to=7914231,by=1989); samples[1]=1

#The gigantic forloop that prints out everything
for(i in 1:length(enrich_CoAnnotation_list)){
  dat=enrich_CoAnnotation_list[[i]]
  dat$xaxis=1:7914231
  dat=dat[samples,]
  dat=melt(dat,id.vars="xaxis")
  
  
  ##Complete graph
  #library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  #ggplot(data=dat,aes(x=xaxis,y=value,color=variable))+
  #  geom_line()+
  #  theme(text = element_text(size=20))+
  #  geom_hline(yintercept=1, linetype=2, color="red", size=1)+
  #  labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text')
  
  ##Graph that contains only combinations of annotation sets
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  dat2=dat[dat$variable %in% c("any","anyTwo","anyThree","anyFour","anyFive"),]
  ggplot(data=dat2,aes(x=xaxis,y=value,color=variable))+
    scale_y_continuous(limits = c(0, 16000))+
    geom_line(size=1,alpha=0.8)+
    theme_minimal()+
    theme(text = element_text(size=20))+
    geom_hline(yintercept=1, linetype=2, color="red", size=1)+
    labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text')+
    guides(color=guide_legend("Annotation Set"))
  ggsave(paste(distances[i],"_combination.png",sep=""),width=12,height=8)
  
  ##Graph that contains different annotation sets
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  dat3=dat[dat$variable %in% c("Pwy","pcomplex","operon","regulator","kegg_modules",
                               "Pwy_and_kegg_modules","Pwy_or_kegg_modules"),]
  ggplot(data=dat3,aes(x=xaxis,y=value,color=variable))+
    scale_y_continuous(limits = c(0, 16000))+
    geom_line(size=1,alpha=0.8)+
    theme_minimal()+
    theme(text = element_text(size=20))+
    geom_hline(yintercept=1, linetype=2, color="red", size=1)+
    labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text')+
    guides(color=guide_legend("Annotation Set"))
  ggsave(paste(distances[i],"_diffAnnotSet.png",sep=""),width=12,height=8)
  
}

