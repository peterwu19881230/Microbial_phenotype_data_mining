#Parallel programming ref:
#https://blog.dominodatalab.com/simple-parallelization/
#https://cran.r-project.org/web/packages/foreach/vignettes/nested.pdf

library(foreach)
library(doParallel)
library(parallel)

numCores=detectCores()-1
cl=makeCluster(numCores)
#Note: For windows I also need this: https://www.r-bloggers.com/parallel-r-loops-for-windows-and-linux/


set.seed(102)
dat=matrix(round(rnorm(3*2)),ncol=3)

#example 1
result=foreach(i=1:dim(dat)[1]) %dopar%{
  sum(dat[i,])
}

#example 2
result2=foreach(i=1:dim(dat)[1]) %:% 
  foreach(j=1:dim(dat)[2]) %dopar% {
    dat[i,j]^2
  }

#example 3
result3=foreach(i=1:dim(dat)[1]) %dopar% {
  library(foreach)
  dat=dat+1
  foreach(j=1:dim(dat)[2]) %dopar% {
    dat[i,j]^2
  }
}

stopCluster(cl)
