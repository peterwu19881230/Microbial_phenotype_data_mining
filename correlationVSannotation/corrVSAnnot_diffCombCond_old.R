#Find the most useful conditions that give the best corrVSannot


#This function takes a distance object as input and output a melted dataframe
#!(NA will be kicked out)
melt_dist=function(distance){
  m1=as.matrix(distance)
  m1[upper.tri(m1)]=NA #I suppose that no distance can be NA, so I can use this to do filtering
  diag(m1)=NA
  library(reshape2)
  m2=melt(m1) #I suppose that no distance can be NA, so I can use this to do filtering
  if(class(m2$Var1)!="character") m2$Var1=as.character(m2$Var1) #This is to prevent numeric names (the class should still be "character") being converted to "numeric" by melt()
  if(class(m2$Var2)!="character") m2$Var2=as.character(m2$Var2) #This is to prevent numeric names (the class should still be "character") being converted to "numeric" by melt()
  
  m2=m2[!is.na(m2[,3]) | is.nan(m2[,3]),] # "| is.nan(m2[,3])" is used because I want to keep the NaN values
  
  m2=m2[,c(2,1,3)] ##reorder the columns
  names(m2)=c("Object 1","Object 2","Distance") ##name the columns
  return(m2)
}


# Exhasutively get all combitations of conditions for no. of conditions= 3 ~ 324
'
set.seed(101)
dat=matrix(rnorm(10*5),ncol=5) #change to All_Data_NAimputed later
n=5 #change to 324 later
annotation=c(rep(1,10),rep(0,35)) #change to any annotation T or F
strain1strain2PCC=list()
for(i in 3:n){
  comb=t(combn(n,i))
  strain1strain2PCC_=list()
  strain1strain2PCC_$comb=comb
  strain1strain2PCC_$dat=t(combn(dim(dat)[1],2))
  for(j in 1:dim(comb)[1]){
    strain1strain2PCC_$dat=cbind(strain1strain2PCC_$dat,melt_dist(1-abs(cor(t(dat[,comb[j,]]))))[,3]) #melt_dist is defined above
  }
  strain1strain2PCC_$dat=cbind(strain1strain2PCC_$dat,annotation)
  strain1strain2PCC[[i-2]]=strain1strain2PCC_
}

# Determine the total no. of co-annotations at cutoff |PCC|=0.6 (PCC based distance = 0.4)
result=lapply(strain1strain2PCC,FUN = function(comb_dat){
  apply(as.matrix(comb_dat$dat[,3:(dim(comb_dat$dat)[2]-1)]),2,FUN=function(pcc_cols){
    sum(ifelse(pcc_cols<=0.4 & comb_dat$dat[,"annotation"]==1,1,0))
  })

})
'
#parallelized version
#Ref: 
#https://cran.r-project.org/web/packages/foreach/vignettes/nested.pdf
#https://www.r-bloggers.com/the-wonders-of-foreach/
library(parallel); library(doParallel); library(foreach); library(doSNOW);

start=Sys.time()

numCores=detectCores()-1
cl = makeCluster(numCores, type = "SOCK")
registerDoSNOW(cl)


# Exhasutively get all combitations of conditions for no. of conditions= 3 ~ 324
set.seed(101)
n=5 #change to 324 later
dat=matrix(rnorm(10*n),ncol=n) #change to All_Data_NAimputed later
annotation=c(rep(1,10),rep(0,35)) #change to any annotation T or F
strain1strain2PCC=list()
for(i in 3:n){
  comb=t(combn(n,i))
  strain1strain2PCC_=list()
  strain1strain2PCC_$comb=comb
  
  strain1strain2PCC_$dat=cbind(t(combn(dim(dat)[1],2)),
                               foreach(j=1:dim(comb)[1], .combine="cbind") %dopar% {
                                 melt_dist(1-abs(cor(t(dat[,comb[j,]]))))[,3] #melt_dist is defined above
                               }
  )
  
  strain1strain2PCC_$dat=cbind(strain1strain2PCC_$dat,annotation)
  strain1strain2PCC[[i-2]]=strain1strain2PCC_
}

# Determine the total no. of co-annotations at cutoff |PCC|=0.6 (PCC based distance = 0.4)
result=lapply(strain1strain2PCC,FUN = function(comb_dat){
  apply(as.matrix(comb_dat$dat[,3:(dim(comb_dat$dat)[2]-1)]),2,FUN=function(pcc_cols){
    sum(ifelse(pcc_cols<=0.4 & comb_dat$dat[,"annotation"]==1,1,0))
  })
  
})

stopCluster(cl)

end=Sys.time()
end-start

#save(strain1strain2PCC,result,file="diffCombCond.RData")



