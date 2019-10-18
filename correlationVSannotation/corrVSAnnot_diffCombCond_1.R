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


#

#parallelized version
#Ref: 
#https://cran.r-project.org/web/packages/foreach/vignettes/nested.pdf
#https://www.r-bloggers.com/the-wonders-of-foreach/
library(parallel); library(doParallel); library(foreach); library(doSNOW);

start=Sys.time()

numCores=detectCores()-1
cl = makeCluster(numCores, type = "SOCK")
registerDoSNOW(cl)




