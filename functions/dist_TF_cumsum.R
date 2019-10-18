# This functions takes a named and tagged dataset (a data frame) and output how well some differnet distances are in evaluating the dataset. 
# melt_dist(), another self-written function, is required. 
# (!This function doesn't deal with many-to-many relationship yet (tag-name) )
# (!This function doesn't tolerate NAs in the data variable)

dist_TF_cumsum=function(data,attribute_list,dist_metric){
  if(class(rownames(data))!="character")rownames(data)=as.character(rownames(data))
  
  
  distance=dist_metric(data)
  melted=melt_dist(distance)
  
  distance_table=melted
  
  
  names(distance_table)[3:(dim(distance_table)[2])]=c("distance")
  
  #===According to profiler, this block is the bottle neck=============================
  sameORnot=c()
  for(i in 1:dim(distance_table)[1]){
    id1=distance_table[i,1]
    id2=distance_table[i,2]
    
    sameORnot[i]=ifelse(  anyNA(attribute_list[[id1]]) || 
                            anyNA(attribute_list[[id2]]) ||
                            sum( attribute_list[[id1]] %in% attribute_list[[id2]])==0,  #The way attribute_list[[id1]] works is like attribute_list$id1, where id1 is a string. I haven't tested the function when id1 is a number
                          0,1)
  }
  
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



#Below was the older algorithm (I replaced these 4 lines with the above 4 lines)
#preliminary testing shows that the above new algorithm is more than 2 times faster than the older version

#sameORnot[i]=ifelse( sum( attribute_list[[id1]] %in% attribute_list[[id2]])>=1 & #The way attribute_list[[id1]] works is like attribute_list$id1, where id1 is a string. I haven't tested the function when id1 is a number 
#                       !(anyNA(attribute_list[[id1]])) & 
#                       !(anyNA(attribute_list[[id2]])),
#                     1,0)




#Beolow is another approach, but it is 2 times slower than the current algorithm
#Ref: https://stackoverflow.com/questions/24614391/intersect-all-possible-combinations-of-list-elements

#replacement code for the code in the bottle neck block==================================================

#ll <- combn( attr_list , 2 , simplify = FALSE )

#TF=sapply( ll , function(x){
#  ifelse( (NA %in% intersect( x[[1]] , x[[2]] )) ||
#            identical(intersect( x[[1]] , x[[2]] ), logical(0)) ||
#            identical(intersect( x[[1]] , x[[2]] ), numeric(0)) ||
#            identical(intersect( x[[1]] , x[[2]] ), character(0)),
#          0,1)
#})

#Obj1_Obj2_TF=data.frame(
#  (purrr::map(ll,~names(.x)[1]) %>% unlist),
#  (purrr::map(ll,~names(.x)[2]) %>% unlist),
#  sameORnot=TF
#)

#names(Obj1_Obj2_TF)=c("Object 1","Object 2","sameORnot")

#distance_table_sameORnot=merge(distance_table,Obj1_Obj2_TF)

#======================================================================



#Toy data for testing:
#attr_list=list(A=c("one", "two", "three", "four"), B=c("one", "two"), C=c("two", "four", "five", "six"), D=c("six", "seven"),E=NA,f=NA,G=NA)
#dat=data.frame(A=rnorm(1000),B=rnorm(1000),C=rnorm(1000),D=rnorm(1000),E=rnorm(1000),f=rnorm(1000),G=rnorm(1000)) %>% t



