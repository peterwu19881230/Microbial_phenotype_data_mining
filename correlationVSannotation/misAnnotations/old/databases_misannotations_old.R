class(id_allAttributes)
id_pwyAnnot=id_allAttributes[,c("ids","Pwy")] %>% unique
dim(id_pwyAnnot)



#The following mis-annotation step might not be the best (but the point is to somehow contaminate the database, and the following might be good enough)
#Note: this way the no. of misannotations is underestimated because for example: sum(c(1,2,NA,4,5)!=c(2,3,4,5,6),na.rm=T) gives 4 while there are actually 5 differences

#Difference after shuffling annotations


#Differnce between the origianl and complete random
fake_pwyAnnot=id_pwyAnnot
sum(!is.na(fake_pwyAnnot$Pwy)) #total annotations= 2344

set.seed(102)
random_index=base::sample(dim(fake_pwyAnnot)[1],dim(fake_pwyAnnot)[1])
sum(fake_pwyAnnot$Pwy!=fake_pwyAnnot$Pwy[random_index],na.rm=T) #975 annotations are different when they are re-ordered randomly
#=> This indicates no matter how many times I shuffle the annotations, it's hard to get more than 975 annotations (using different seed gives me a no. at around 1000)

start=Sys.time()

fake_pwyAnnot=id_pwyAnnot
index=1:dim(id_pwyAnnot)[1]


i=1
j=1
repeat{
  set.seed(i)
  toBeSwitched=base::sample(dim(id_pwyAnnot)[1],2)
  index[toBeSwitched[1]]=toBeSwitched[2];index[toBeSwitched[2]]=toBeSwitched[1]
  
  if(sum(fake_pwyAnnot$Pwy!=fake_pwyAnnot$Pwy[index],na.rm=T)>=100*j){ #maybe change to "==100*j" later?
    fake_pwyAnnot[,j+2]=fake_pwyAnnot$Pwy[index]
    print(paste("no. ",j," mis-annotation set is done",sep=""))
    j=j+1
  } 
  
  if(sum(fake_pwyAnnot$Pwy!=fake_pwyAnnot$Pwy[index],na.rm=T)==1000){ 
    ##Why 1000?: no matter how many times I shuffle the annotations, it's hard to get more than 975 annotations (using different seed gives me a no. at around 1000)
    break
  }
  
  i=i+1
}

end=Sys.time()
end-start #Time difference of 6.437369 secs





#plot to show how different columns contain various no. of mis-annotations
if(!file.exists("Data/mis_annotation_dfs_old.RData")){  #If the file cannot be found, re-process
  
  no_misAnnot=c()
  for(i in 3:dim(fake_pwyAnnot)[2]){
    no_misAnnot[i-2]=sum(fake_pwyAnnot$Pwy!=fake_pwyAnnot[,i],na.rm=T)
  }
  no_misAnnot #101  200  301  400  501  600  700  800  900 1000
  
  
  dat=All_Data_NAimputed
  
  
  start.time = Sys.time()
  
  mis_annotation_dfs_old=list()
  for(i in 2:dim(fake_pwyAnnot)[2]){
    attribute_list=attr_list(fake_pwyAnnot$ids,fake_pwyAnnot[,i])
    mis_annotation_dfs_old[[i-1]]=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = pcc_dist)
  }
  save(mis_annotation_dfs_old,file="Data/mis_annotation_dfs_old.RData")
  
  
  
  
  end.time = Sys.time()
  end.time - start.time #Time difference of 3.164703 hours (for my PC: Time difference of 1.985411 hours)
  
}else load("Data/mis_annotation_dfs_old.RData")






#pcc at no. 1000th, 2000th, 3000th
result=mis_annotation_dfs_old

1-result[[1]][[1]]$distance[c(1000,2000,3000)] #0.7886812 0.6652216 0.5774545


#Plot the result:

##1. Line plot
samples=1:4000
#samples=seq(from=0,to=7914231,by=1989); samples[1]=1

library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+theme(text = element_text(size=20))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[1]][[1]]$cumsum[samples]/((sum(result[[1]][[1]]$sameORnot)/length(result[[1]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="0"))+ #0=original (use 0 instead of "original to prevent wrong re-ordering of the side bar")
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[2]][[1]]$cumsum[samples]/((sum(result[[2]][[1]]$sameORnot)/length(result[[2]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="1"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[3]][[1]]$cumsum[samples]/((sum(result[[3]][[1]]$sameORnot)/length(result[[3]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="2"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[4]][[1]]$cumsum[samples]/((sum(result[[4]][[1]]$sameORnot)/length(result[[4]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="3"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[5]][[1]]$cumsum[samples]/((sum(result[[5]][[1]]$sameORnot)/length(result[[5]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="4"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[6]][[1]]$cumsum[samples]/((sum(result[[6]][[1]]$sameORnot)/length(result[[6]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="5"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[7]][[1]]$cumsum[samples]/((sum(result[[7]][[1]]$sameORnot)/length(result[[7]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="6"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[8]][[1]]$cumsum[samples]/((sum(result[[8]][[1]]$sameORnot)/length(result[[8]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="7"))+
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[9]][[1]]$cumsum[samples]/((sum(result[[9]][[1]]$sameORnot)/length(result[[9]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="8"))+  
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[10]][[1]]$cumsum[samples]/((sum(result[[10]][[1]]$sameORnot)/length(result[[10]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="9"))+  
  geom_line(data=data.frame(no=samples,foldEnrichment=result[[11]][[1]]$cumsum[samples]/((sum(result[[11]][[1]]$sameORnot)/length(result[[11]][[1]]$sameORnot))*samples)),aes(no,foldEnrichment,color="ten"))+ #ten is to prevent wrong re-ordering of the side bar
  geom_hline(yintercept=1, linetype=2, color="red", size=1)+
  labs(x = "Ranking of Distance",y="Fold Enrichment (normalized by random expectation)",aesthetic='custom text')



##2. Corresponding barplots (at rank 10,20,50,100,200,500,1000,2000,3000)
samples=c(10,20,50,100,200,500,1000,2000,3000)
dat=data.frame(no=samples,
               foldEnrichment=result[[1]][[1]]$cumsum[samples]/((sum(result[[1]][[1]]$sameORnot)/length(result[[1]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[2]][[1]]$cumsum[samples]/((sum(result[[2]][[1]]$sameORnot)/length(result[[2]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[3]][[1]]$cumsum[samples]/((sum(result[[3]][[1]]$sameORnot)/length(result[[3]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[4]][[1]]$cumsum[samples]/((sum(result[[4]][[1]]$sameORnot)/length(result[[4]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[5]][[1]]$cumsum[samples]/((sum(result[[5]][[1]]$sameORnot)/length(result[[5]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[6]][[1]]$cumsum[samples]/((sum(result[[6]][[1]]$sameORnot)/length(result[[6]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[7]][[1]]$cumsum[samples]/((sum(result[[7]][[1]]$sameORnot)/length(result[[7]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[8]][[1]]$cumsum[samples]/((sum(result[[8]][[1]]$sameORnot)/length(result[[8]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[9]][[1]]$cumsum[samples]/((sum(result[[9]][[1]]$sameORnot)/length(result[[9]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[10]][[1]]$cumsum[samples]/((sum(result[[10]][[1]]$sameORnot)/length(result[[10]][[1]]$sameORnot))*samples),
               foldEnrichment=result[[11]][[1]]$cumsum[samples]/((sum(result[[11]][[1]]$sameORnot)/length(result[[11]][[1]]$sameORnot))*samples)
               )


#10
p=list()
for(i in 1:9){
  dat_i=data.frame(noOfMisAnnot=c("original",
                                   "101",
                                   "200",
                                   "301",
                                   "400",
                                   "501",
                                   "600",
                                   "700",
                                   "800",
                                   "900",
                                   "1000"
  ),foldEnrichment=as.vector(t(dat[i,2:12])))
  
  dat_i$noOfMisAnnot=factor(dat_i$noOfMisAnnot,levels=dat_i$noOfMisAnnot) #change x-axis to a factor variable to prevent ggplot's automatic re-ordering
  
  p[[i]]=ggplot(data=dat_i,aes(x=noOfMisAnnot,y=foldEnrichment))+geom_bar(stat="identity")+
          labs(title=paste("No. of ranking = ",samples[i],sep=""),x="No. of mis-annotations",y="Fold enrichment")
           
}




gridExtra::grid.arrange(grobs=p,ncol=3)




