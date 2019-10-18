#Parallel programming tested on ADA
#ssh peterwu19881230@ada.tamu.edu
#module load R/3.5.0-iomkl-2017b-recommended-mt


##Somthing below cause error on ADA (and will run for a very long time)
##if(!require(pacman)){
##  install.packages("pacman")
##  library(pacman)
##}
##pacman::p_load(tidyverse,xlsx,factoextra,pheatmap,ComplexHeatmap) #Note ggplot2 is loaded when factoextra is loaded 

#Test 1
set.seed(0)
largeData=matrix(rnorm(10^4*10^2*10^2),nrow=10^2)

start.time=Sys.time()
nonpar=apply(largeData,MARGIN=1,FUN=median) 
end.time=Sys.time()
end.time-start.time



#Tutorial on parallel apply family: http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
library(parallel)
no_cores <- detectCores()-1 
cl <- makeCluster(no_cores)
start.time=Sys.time()
par=parApply(cl=cl,X=largeData,MARGIN=1,FUN=median)
end.time=Sys.time()
end.time-start.time


#Test 2
set.seed(0)
largeData=rnorm(7914231)


n=c(3:22,27,28,31,41,47,48)

start.time = Sys.time()
set.seed(101)
randomPCC=list()
i=1
for(num in n){
  randomPCC[[i]]=sapply(1:5000,FUN=function(iteration){
    mean(abs(sample(largeData,num)))
  })
  i=i+1
}
end.time = Sys.time()
end.time - start.time #Time difference of 9.82577 mins (When I run R command line on ADA)




library(parallel)
no_cores <- detectCores()-1 
cl <- makeCluster(no_cores)

##Must use clusterExort or specify within the anonymous function: 
##Ref: https://stackoverflow.com/questions/10095956/parsapply-not-finding-objects-in-global-environment
clusterExport(cl,"largeData") #If multiple objs need to be recognized, concatnate the input to be: c("Obj1","Obj2","Obj3"...)

start.time = Sys.time()
set.seed(101)
randomPCC=list()
i=1
for(num in n){
  clusterExport(cl,"num")
  randomPCC[[i]]=parSapply(cl=cl,X=1:5000,FUN=function(iteration){
    mean(abs(sample(largeData,num)))
  })
  i=i+1
}
end.time = Sys.time()
end.time - start.time #Time difference of 1.44726 mins (When I run R command line on ADA)




