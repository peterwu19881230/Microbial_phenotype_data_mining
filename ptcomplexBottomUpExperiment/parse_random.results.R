#Goal: Parse the results of the random experiment

## The No. of genes in the pathways used
## j=2 is to deal with situation where 3 or more genes are in the pathway but only 2 genes are found in Nichols' 
j=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,26,27,28,31,41,48)


## Read each processed data and combind them into 1 dataframe
## (!) For larger number of j, the random fraction isn't just 1/number. It would be large (eg. 1/48(~0.02) is not equal to 0.03815026, the calculated fraction )
pwyBottomUpExp.random=data.frame()
for(i in j){
  file.name=paste("avg.random_tables.j",i,".RData",sep="")
  path=paste("pwyBottomToTopExperiment/random.result.FromADA/",file.name,sep="")
  load(path)
  pwyBottomUpExp.random=rbind(pwyBottomUpExp.random, avg.random_tables)
  
}

save(pwyBottomUpExp.random,file="Data/pwyBottomUpExp.random.RData")