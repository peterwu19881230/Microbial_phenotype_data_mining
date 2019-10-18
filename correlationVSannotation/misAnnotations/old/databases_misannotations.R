#Assess the effect of misannotations on correlation VS. annotation experiment
##1. Shuffled annotations
##2. With some annotations thrown out
##3. compare results of different quantitative distances 
##4. compare results of different qualitative distances


##1. Shuffled annotations
load("Data/mis_annotation_random.RData")
load("Data/mis_annotation_quantitative_list.RData")


#pcc at no. 1000th, 2000th, 3000th
result=mis_annotation_quantitative_list
random_df=mis_annotation_random

1-result[[1]][[1]]$distance[c(1000,2000,3000)] #0.7886812 0.6652216 0.5774545



#Plot the result:

##1. Line plot
samples=1:4000
#samples=1:100
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1

library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+theme(text = element_text(size=20))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[1]][[1]][[1]]$cumsum[samples]/((sum(result[[1]][[1]][[1]]$sameORnot)/length(result[[1]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="0"))+ #0=original (use 0 instead of "original to prevent wrong re-ordering of the side bar")
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[4]][[1]][[1]]$cumsum[samples]/((sum(result[[4]][[1]][[1]]$sameORnot)/length(result[[4]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="100"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[6]][[1]][[1]]$cumsum[samples]/((sum(result[[6]][[1]][[1]]$sameORnot)/length(result[[6]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="1000"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=random_df[[1]]$cumsum[samples]/((sum(random_df[[1]]$sameORnot)/length(random_df[[1]]$sameORnot))*samples)),aes(no,foldEnrichment),linetype=2, color="red", size=1)+ 
  labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text')





##2. Corresponding barplots (at rank 10,20,50,100,200,500,1000,2000,3000)
samples=c(10,20,50,100,200,500,1000,2000,3000)
dat=data.frame(no=samples,
               foldEnrichment=result[[1]][[1]][[1]]$cumsum[samples]/((sum(result[[1]][[1]][[1]]$sameORnot)/length(result[[1]][[1]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[4]][[1]][[1]]$cumsum[samples]/((sum(result[[4]][[1]][[1]]$sameORnot)/length(result[[4]][[1]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[6]][[1]][[1]]$cumsum[samples]/((sum(result[[6]][[1]][[1]]$sameORnot)/length(result[[6]][[1]][[1]]$sameORnot))*samples),
               foldEnrichment=random_df[[1]]$cumsum[samples]/((sum(random_df[[1]]$sameORnot)/length(random_df[[1]]$sameORnot))*samples)
)


#10
p=list()
for(i in seq_along(samples)){
  dat_i=data.frame(noOfMisAnnot=c("original",
                                  "100",
                                  "1000",
                                  "random"
  ),foldEnrichment=as.vector(t(dat[i,2:dim(dat)[2]])))
  
  dat_i$noOfMisAnnot=factor(dat_i$noOfMisAnnot,levels=dat_i$noOfMisAnnot) #change x-axis to a factor variable to prevent ggplot's automatic re-ordering
  
  p[[i]]=ggplot(data=dat_i,aes(x=noOfMisAnnot,y=foldEnrichment))+geom_bar(stat="identity")+
    labs(title=paste("No. of ranking = ",samples[i],sep=""),x="No. of mis-annotations",y="Fold enrichment")
  
}




gridExtra::grid.arrange(grobs=p,ncol=3)



##2. With some annotations thrown out

load("Data/mis_annotation_random.RData")
load("Data/thrown_annotation_quantitative_list.RData")


#pcc at no. 1000th, 2000th, 3000th
result=thrown_annotation_quantitative_list
random_df=mis_annotation_random

1-result[[1]][[1]]$distance[c(1000,2000,3000)] #0.7886812 0.6652216 0.5774545



#Plot the result:

##1. Line plot
samples=1:4000
#samples=1:100
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1

library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+theme(text = element_text(size=20))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[1]][[1]][[1]]$cumsum[samples]/((sum(result[[1]][[1]][[1]]$sameORnot)/length(result[[1]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="0"))+ #0=original (use 0 instead of "original to prevent wrong re-ordering of the side bar")
  #geom_line(data=data.frame(no=samples,foldEnrichment=result[[2]][[1]][[1]]$cumsum[samples]/((sum(result[[2]][[1]][[1]]$sameORnot)/length(result[[2]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="1"))+
  #geom_line(data=data.frame(no=samples,foldEnrichment=result[[3]][[1]][[1]]$cumsum[samples]/((sum(result[[3]][[1]][[1]]$sameORnot)/length(result[[3]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="10"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[4]][[1]][[1]]$cumsum[samples]/((sum(result[[4]][[1]][[1]]$sameORnot)/length(result[[4]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="100"))+
  #geom_line(data=data.frame(no=samples,foldEnrichment=result[[5]][[1]][[1]]$cumsum[samples]/((sum(result[[5]][[1]][[1]]$sameORnot)/length(result[[5]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="500"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[6]][[1]][[1]]$cumsum[samples]/((sum(result[[6]][[1]][[1]]$sameORnot)/length(result[[6]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="1000"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[7]][[1]][[1]]$cumsum[samples]/((sum(result[[7]][[1]][[1]]$sameORnot)/length(result[[7]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="1500"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[8]][[1]][[1]]$cumsum[samples]/((sum(result[[8]][[1]][[1]]$sameORnot)/length(result[[8]][[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="2000"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=random_df[[1]]$cumsum[samples]/((sum(random_df[[1]]$sameORnot)/length(random_df[[1]]$sameORnot))*samples)),aes(no,foldEnrichment),linetype=2, color="red", size=1)+ 
  labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text')




#old way of plotting
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+theme(text = element_text(size=20))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[1]][[1]][[1]]$cumsum[samples]),aes(no,foldEnrichment,color="0"))+ #0=original (use 0 instead of "original to prevent wrong re-ordering of the side bar")
  #geom_line(data=data.frame(no=samples,foldEnrichment=result[[2]][[1]][[1]]$cumsum[samples]),aes(no,foldEnrichment,color="1"))+
  #geom_line(data=data.frame(no=samples,foldEnrichment=result[[3]][[1]][[1]]$cumsum[samples]),aes(no,foldEnrichment,color="10"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[4]][[1]][[1]]$cumsum[samples]),aes(no,foldEnrichment,color="100"))+
  #geom_line(data=data.frame(no=samples,foldEnrichment=result[[5]][[1]][[1]]$cumsum[samples]),aes(no,foldEnrichment,color="500"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[6]][[1]][[1]]$cumsum[samples]),aes(no,foldEnrichment,color="1000"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[7]][[1]][[1]]$cumsum[samples]),aes(no,foldEnrichment,color="1500"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[8]][[1]][[1]]$cumsum[samples]),aes(no,foldEnrichment,color="2000"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=random_df[[1]]$cumsum[samples]/((sum(random_df[[1]]$sameORnot)/length(random_df[[1]]$sameORnot))*samples)),aes(no,foldEnrichment),linetype=2, color="red", size=1)+ 
  labs(x = "Ranking of Distance",y="Cumulative sum",aesthetic='custom text')



#additonal plots to see whether "fold enrichment" should be a good y axis

result=thrown_annotation_quantitative_list_2
random_df=mis_annotation_random

#Plot the result:

##1. Line plot
samples=1:200
#samples=1:100
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1


n=c(1,1:20*100)

foldEnrichment.df=data.frame()
for(i in seq(thrown_annotation_quantitative_list_2)){ 
  foldEnrichment.df=rbind(
    foldEnrichment.df,
    data.frame(no=samples,
               foldEnrichment=result[[i]][[1]][[1]]$cumsum[samples]/((sum(result[[i]][[1]][[1]]$sameORnot)/length(result[[i]][[1]][[1]]$sameORnot))*samples),
               no_of_thrownouts=as.character(n[i]))
  )
}


#show all the thrown-out annotations
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
color_function=colorRampPalette(c("blue","green","red","pink"))

ggplot()+theme_minimal()+theme(text = element_text(size=20))+
  geom_line(data=foldEnrichment.df,aes(x=no,y=foldEnrichment,color=no_of_thrownouts),alpha=0.8,size=1)+ 
  geom_line(data=data.frame(no=samples,foldEnrichment=random_df[[1]]$cumsum[samples]/((sum(random_df[[1]]$sameORnot)/length(random_df[[1]]$sameORnot))*samples)),aes(no,foldEnrichment),linetype=2, color="black", size=1)+ 
  labs(x = "Ranking of Distance (PCC)",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text')+
  scale_color_manual(values=color_function(21))









##3. compare results of different quantitative distances




##1. Line plot
samples=1:4000
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1

#Using original EcoCyc pwy annotation set
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme(text = element_text(size=20))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result$pearson[[1]]$cumsum[samples]/((sum(result$pearson[[1]]$sameORnot)/length(result$pearson[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Pearson"))+ 
  geom_line(data=data.frame(no=samples,foldEnrichment=result$spearman[[1]]$cumsum[samples]/((sum(result$spearman[[1]]$sameORnot)/length(result$spearman[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Spearman"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result$euclidean[[1]]$cumsum[samples]/((sum(result$euclidean[[1]]$sameORnot)/length(result$euclidean[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Euclidean"))+
  geom_hline(yintercept=1, linetype=2, color="red", size=1)+
  labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text',color="Distance")





#Using 500 mis-annotatinos
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme(text = element_text(size=20))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result_500$pearson[[1]]$cumsum[samples]/((sum(result_500$pearson[[1]]$sameORnot)/length(result_500$pearson[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Pearson"))+ 
  geom_line(data=data.frame(no=samples,foldEnrichment=result_500$spearman[[1]]$cumsum[samples]/((sum(result_500$spearman[[1]]$sameORnot)/length(result_500$spearman[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Spearman"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result_500$euclidean[[1]]$cumsum[samples]/((sum(result_500$euclidean[[1]]$sameORnot)/length(result_500$euclidean[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Euclidean"))+
  geom_hline(yintercept=1, linetype="dashed", color="red", size=1)+
  labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text',color="Distance")



##combined plot
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+
  theme(text = element_text(size=20))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result$pearson[[1]]$cumsum[samples]/((sum(result$pearson[[1]]$sameORnot)/length(result$pearson[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Pearson"))+ 
  geom_line(data=data.frame(no=samples,foldEnrichment=result$spearman[[1]]$cumsum[samples]/((sum(result$spearman[[1]]$sameORnot)/length(result$spearman[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Spearman"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result$euclidean[[1]]$cumsum[samples]/((sum(result$euclidean[[1]]$sameORnot)/length(result$euclidean[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Euclidean"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result_500$pearson[[1]]$cumsum[samples]/((sum(result_500$pearson[[1]]$sameORnot)/length(result_500$pearson[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Pearson"),linetype = "dashed")+ 
  geom_line(data=data.frame(no=samples,foldEnrichment=result_500$spearman[[1]]$cumsum[samples]/((sum(result_500$spearman[[1]]$sameORnot)/length(result_500$spearman[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Spearman"),linetype = "dashed")+
  geom_line(data=data.frame(no=samples,foldEnrichment=result_500$euclidean[[1]]$cumsum[samples]/((sum(result_500$euclidean[[1]]$sameORnot)/length(result_500$euclidean[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Euclidean"),linetype = "dashed")+
  geom_hline(yintercept=1, linetype=2, color="red", size=1)+
  labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text',color="Distance")


##(1000 mis-annotation) combined plot
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+
  theme(text = element_text(size=20))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result$pearson[[1]]$cumsum[samples]/((sum(result$pearson[[1]]$sameORnot)/length(result$pearson[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Pearson"))+ 
  geom_line(data=data.frame(no=samples,foldEnrichment=result$spearman[[1]]$cumsum[samples]/((sum(result$spearman[[1]]$sameORnot)/length(result$spearman[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Spearman"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result$euclidean[[1]]$cumsum[samples]/((sum(result$euclidean[[1]]$sameORnot)/length(result$euclidean[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Euclidean"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result_1000$pearson[[1]]$cumsum[samples]/((sum(result_1000$pearson[[1]]$sameORnot)/length(result_1000$pearson[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Pearson"),linetype = "dashed")+ 
  geom_line(data=data.frame(no=samples,foldEnrichment=result_1000$spearman[[1]]$cumsum[samples]/((sum(result_1000$spearman[[1]]$sameORnot)/length(result_1000$spearman[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Spearman"),linetype = "dashed")+
  geom_line(data=data.frame(no=samples,foldEnrichment=result_1000$euclidean[[1]]$cumsum[samples]/((sum(result_1000$euclidean[[1]]$sameORnot)/length(result_1000$euclidean[[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="Euclidean"),linetype = "dashed")+
  geom_hline(yintercept=1, linetype=2, color="red", size=1)+
  labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text',color="Distance")


##4. compare results of different qualitative distances

##After making sure the above graphs will be used I will finish this part













