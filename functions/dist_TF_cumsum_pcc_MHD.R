
#dependent on TernaryConvert(), modified_hamming_distance_3(), melt_dist()
dist_TF_cumsum_pcc_MHD=function(data,ternary=F,thresh=NULL,attribute_list){
  
  if(class(rownames(data))!="character")rownames(data)=as.character(rownames(data))
  
  if(ternary){
    pairwise_cor=cor(t(data)) #Have to transpose to get correlations of the rows
    pcc_dist=1-abs(pairwise_cor) #I use abs() because anti-correlation is better than not having any relationships
    pcc_melted=melt_dist(pcc_dist)
    
    
    ternaryData=TernaryConvert(data,thresh = thresh)
    MHD3_dist=modified_hamming_distance_3(ternaryData)
    MHD3_melted=melt_dist(MHD3_dist)
    
    distance_table=cbind(pcc_melted,MHD3_melted$Distance)
    
    
    names(distance_table)[3:(dim(distance_table)[2])]=c("pcc","MHD3")
    
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
  }else{
    pairwise_cor=cor(t(data)) #Have to transpose to get correlations of the rows
    pcc_dist=1-abs(pairwise_cor) #I use abs() because anti-correlation is better than not having any relationships
    pcc_melted=melt_dist(pcc_dist)
    
    
    distance_table=pcc_melted
    
    
    names(distance_table)[3:(dim(distance_table)[2])]=c("pcc")
    
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
  
}

