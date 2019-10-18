# This functions takes a named and tagged dataset and output how well some differnet distances are in evaluating the dataset. 
# melt_dist(), another self-written function, is required. Here I put it at the bottom
# (!This function doesn't deal with many-to-many relationship yet (tag-name) )
# (!This function doesn't tolerate NAs in the data variable)
dist_TF_cumsum_old=function(data,attribute_list){
  
  if(class(rownames(data))!="character")rownames(data)=as.character(rownames(data))
  
  pairwise_cor=cor(t(data)) #Have to transpose to get correlations of the rows
  pcc_dist=1-abs(pairwise_cor) #I use abs() because anti-correlation is better than not having any relationships
  pcc_squared_dist=1-pairwise_cor^2
  
  if(!require(factoextra)){
    install.packages("factoextra")
    library(factoextra)
  }
  
  euclidean=get_dist(data,method="euclidean") ##get_dist is an enhanced version of dist()
  maximum=get_dist(data,method="maximum")
  manhattan=get_dist(data,method="manhattan")
  canberra=get_dist(data,method="canberra")
  binary=get_dist(data,method="binary")
  minkowski=get_dist(data,method="minkowski",p=3)
  spearman=1-abs(cor(t(data),method="spearman"))
  #dist.All_Data_NAimputed.kendall=get_dist(All_Data_NAimputed,method="kendall") => This works with smaller data frame. But it takes long here so I haven't generated it
  
  pcc_melted=melt_dist(pcc_dist)
  pcc_squared_melted=melt_dist(pcc_squared_dist)$Distance
  euclidean_melted=melt_dist(euclidean)$Distance
  maximum_melted=melt_dist(maximum)$Distance
  manhattan_melted=melt_dist(manhattan)$Distance
  canberra_melted=melt_dist(canberra)$Distance
  binary_melted=melt_dist(binary)$Distance
  minkowski_melted=melt_dist(minkowski)$Distance
  spearman_melted=melt_dist(spearman)$Distance
  
  
  distance_table=cbind(pcc_melted,
                       pcc_squared_melted,euclidean_melted,maximum_melted,manhattan_melted,canberra_melted,binary_melted,minkowski_melted,spearman_melted)
  
  names(distance_table)[3:(dim(distance_table)[2])]=c("pcc","pcc_squared","euclidean","maximum","manhattan","canberra","binary","minkowski","spearman")
  
  sameORnot=c()
  for(i in 1:dim(distance_table)[1]){
    id1=distance_table[i,1]
    id2=distance_table[i,2]
    
    sameORnot[i]=ifelse( sum( attribute_list[[id1]] %in% attribute_list[[id2]])>=1 & #The way attribute_list[[id1]] works is like attribute_list$id1, where id1 is a string. I haven't tested the function when id1 is a number 
                           !(anyNA(attribute_list[[id1]])) & 
                           !(anyNA(attribute_list[[id2]])),
                         1,0)
  }
  
  distance_table_sameORnot=cbind(distance_table,sameORnot)
  
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









