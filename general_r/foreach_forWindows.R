
#Ref: https://www.r-bloggers.com/parallel-r-loops-for-windows-and-linux/

library(foreach)
library(doParallel)
library(parallel)
library(doSNOW) 



numCores=detectCores()-1
cl=makeCluster(numCores)
registerDoSNOW(cl) 


dat=matrix(rnorm(10*7),ncol=7)
iteration=10e4

start.time=Sys.time()

result=foreach(i = 1:7) %dopar% {
  print(i)
  
  for(j in 1:iteration){
    var(dat[,i])
  }
}

stopCluster(cl)

end.time=Sys.time()
end.time-start.time 



start.time=Sys.time()

for(i in 1:7){
  print(i)
  
  for(j in 1:iteration){
    var(dat[,i])
  }
  
}

end.time=Sys.time()
end.time-start.time 
