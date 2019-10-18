# This functions takes a named and tagged dataset (a data frame) and output how well some differnet distances are in evaluating the dataset. 
# melt_dist(), another self-written function, is required. Here I put it at the bottom
# (!This function doesn't deal with many-to-many relationship yet (tag-name) )
# (!This function doesn't tolerate NAs in the data variable)

dist_TF_cumsum_parallel=function(data,attribute_list,dist_metric){
  if(class(rownames(data))!="character")rownames(data)=as.character(rownames(data))
  
  
  distance=dist_metric(data)
  melted=melt_dist(distance)
  
  distance_table=melted
  
  
  names(distance_table)[3:(dim(distance_table)[2])]=c("distance")
  
  #===According to profiler, this block is the bottle neck=============================
  
  library(parallel)
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  clusterExport(cl=cl, "attribute_list") # I tried to pass more than 1 variable and it worked
  
  
  sameORnot=parApply(cl=cl,X=distance_table,MARGIN=1,FUN=function(distance_table_row){
    id1=distance_table_row[1] #have verified that when subsetting by apply(), this is a vector rather than a 1 row dataframe
    id2=distance_table_row[2]
    
    ifelse(  anyNA(attribute_list[[id1]]) || 
             anyNA(attribute_list[[id2]]) ||
             sum( attribute_list[[id1]] %in% attribute_list[[id2]])==0,  #The way attribute_list[[id1]] works is like attribute_list$id1, where id1 is a string. I haven't tested the function when id1 is a number
             0,1)
  }
  )
  
  stopCluster(cl)
  
  
  distance_table_sameORnot=cbind(distance_table,sameORnot)
  #====================================================================================
  
  result_list=list()
  count=1
  for(i in 3:(dim(distance_table_sameORnot)[2]-1)){
    result=distance_table_sameORnot[,c(1,2,i,dim(distance_table_sameORnot)[2])]
    result_list[[count]]=result[order(result[,3]),]
    count=count+1
  }
  
  result_list=lapply(result_list,FUN=function(table){
    table$cumsum=cumsum(table[,4])
    return(table)
  })
  
  return(result_list)
}





#Toy data for testing:
#attr_list=list(A=c("one", "two", "three", "four"), B=c("one", "two"), C=c("two", "four", "five", "six"), D=c("six", "seven"),E=NA,f=NA,G=NA)
#dat=data.frame(A=rnorm(1000),B=rnorm(1000),C=rnorm(1000),D=rnorm(1000),E=rnorm(1000),f=rnorm(1000),G=rnorm(1000)) %>% t

##Does it give the same result as dist_TF_cumsum() give?
##original=dist_TF_cumsum(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
##paral=dist_TF_cumsum_parallel(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
##identical(original,paral) #TRUE


##Show How faster the parallel version is
'
dat=All_Data_NAimputed
id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique()
attribute_list=attr_list(id_pwy$ids,id_pwy$Pwy)

start.time = Sys.time()
original=dist_TF_cumsum(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
end.time = Sys.time()
end.time - start.time #Time difference of 8.388408 mins

start.time = Sys.time()
paral=dist_TF_cumsum_parallel(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
end.time = Sys.time()
end.time - start.time #Time difference of 4.698854 mins

identical(original,paral) #TRUE
'







