#This function uses MHD3
dist_TF_cumsum_MHD3=function(data,ternary=F,thresh=NULL,attribute_list){
  if(class(rownames(data))!="character") rownames(data)=as.character(rownames(data)) 
  
  
  if(ternary==F) {ternaryData=data #If input is already ternary, just use it
  }else ternaryData=TernaryConvert(data,thresh = thresh) 
  
  
  MHD_mat=modified_hamming_distance_3(ternaryData) 
  MHD_melted=melt_dist(MHD_mat)
  
  
  distance_table=MHD_melted
  
  
  names(distance_table)[3:(dim(distance_table)[2])]=c("MHD3")
  
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

