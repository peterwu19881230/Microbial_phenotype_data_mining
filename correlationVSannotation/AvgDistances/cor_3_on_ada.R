#ssh peterwu19881230@ada.tamu.edu
#module load R/3.5.0-iomkl-2017b-recommended-mt


#Copy and paste the following on ADA:
load("cor_strains.RData")

unique_cor_strains=as.numeric(as.dist(cor_strains))
abs_unique_cor_strains=abs(unique_cor_strains)



library(parallel)
no_cores <- detectCores()-1 
cl <- makeCluster(no_cores)

##Must use clusterExort or specify within the anonymous function: 
##Ref: https://stackoverflow.com/questions/10095956/parsapply-not-finding-objects-in-global-environment


start.time = Sys.time()
#set.seed(101)

clusterExport(cl,"unique_cor_strains")
randomPCC_3=parSapply(cl=cl,X=1:(5000*100),FUN=function(iteration){
    dat=sample(unique_cor_strains,3)
    mean=mean(dat)
    sd=sd(dat)
    return(c(mean,sd))
  })


clusterExport(cl,"abs_unique_cor_strains")
random_absPCC_3=parSapply(cl=cl,X=1:(5000*100),FUN=function(iteration){
  dat=sample(abs_unique_cor_strains,3)
  mean=mean(dat)
  sd=sd(dat)
  return(c(mean,sd))
})


end.time = Sys.time()
end.time - start.time 

save(randomPCC_3,random_absPCC_3,file = "large_cor_3_sampling.RData")

#(from local machine:)
#scp peterwu19881230@ada.tamu.edu:/home/peterwu19881230/large_cor_3_sampling.RData /Users/peterwu/Dropbox/Nichols_Data_mining/Data/sourced



