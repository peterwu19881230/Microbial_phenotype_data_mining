#Parallel processing

set.seed(0)
largeData=matrix(rnorm(10^2*10^2),nrow=10^2)

library(dplyr)
start.time=Sys.time()
nonpar=apply(largeData,MARGIN=1,FUN=median) %>% replicate(n=10^8)
end.time=Sys.time()
end.time-start.time




#Parallel: Much longer than "Not parallel"
#Tutorial on parallel apply family: http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
library(parallel)
no_cores <- detectCores()-2 
cl <- makeCluster(no_cores)
start.time=Sys.time()
par=parApply(cl=cl,X=largeData,MARGIN=1,FUN=median) %>% replicate(n=10^8) 
stopCluster(cl)
end.time=Sys.time()
end.time-start.time

#Conclusion: maybe try on ADA?
