#attribute_list=original_attribute_list[strain_to_merge$ids]

#Method 1: https://stackoverflow.com/questions/24614391/intersect-all-possible-combinations-of-list-elements
l=attribute_list
nms <- combn( names(l) , 2 , FUN = paste0 , collapse = "" , simplify = FALSE )

# Make the combinations of list elements
ll <- combn( l , 2 , simplify = FALSE )

# Intersect the list elements
out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )

# Output with names
setNames( out , nms )



#Method 2

dist_TF_cumsum_v2=function(data,id_annotation,dist_metric){
  if(class(rownames(data))!="character")rownames(data)=as.character(rownames(data))
  
  
  distance=dist_metric(data)
  melted=melt_dist(distance)
  
  distance_table=melted
  
  
  names(distance_table)[3:(dim(distance_table)[2])]=c("distance")
  
  
#================================================================================ I don't understand why this implementation is the slowest
  if(class(id_annotation)!="data.frame"|class(id_annotation[,1])!="character"|class(id_annotation[,2])!="character") stop("id_annotation is not a data frame of 2 character columns")
  
  names(id_annotation)=c("id","annotation")
  annot_moreThanTwo=table(id_annotation$annotation)[table(id_annotation$annotation)>1] %>% names 
  
  coAnnotTable=lapply(annot_moreThanTwo,FUN=function(annot){
    id_annotation$id[id_annotation$annotation==annot] %>% combn(2) %>% t
  }) %>%  BiocGenerics::Reduce(f=rbind) %>% as.data.frame(stringsAsFactors=F) %>% unique
  
  if(dim(coAnnotTable)[1]==0) stop("No co-annotated pairs found")
  
  names(coAnnotTable)=names(distance_table)[1:2]
  coAnnotTable$sameORnot=1
  
  
  distance_table_sameORnot=merge(distance_table,coAnnotTable,all.x=T,all.y=T)
  distance_table_sameORnot$sameORnot[is.na(distance_table_sameORnot$sameORnot)]=0
  
  #================================================================================
  
  
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


#Method 3
dist_TF_cumsum_v3=function(data,attribute_list,dist_metric){
  if(class(rownames(data))!="character")rownames(data)=as.character(rownames(data))
  
  
  distance=dist_metric(data)
  melted=melt_dist(distance)
  
  distance_table=melted
  
  
  names(distance_table)[3:(dim(distance_table)[2])]=c("distance")
  
  sameORnot=c()
  for(i in 1:dim(distance_table)[1]){
    id1=distance_table[i,1]
    id2=distance_table[i,2]
    
    if(anyNA(attribute_list[[id1]]) | anyNA(attribute_list[[id2]])){
      sameORnot[i]=0
      next
    } 
    
    sameORnot[i]=ifelse( sum( attribute_list[[id1]] %in% attribute_list[[id2]])==0,  #The way attribute_list[[id1]] works is like attribute_list$id1, where id1 is a string. I haven't tested the function when id1 is a number 
                         0,1)
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

#Test cases
start.time = Sys.time()
number=500
temp3=dist_TF_cumsum(data=All_Data_NAimputed[1:number,],attribute_list =attr_list(ECK.Pathway$ids,ECK.Pathway$Data) ,dist_metric = pcc_dist)
end.time = Sys.time()
end.time - start.time 

start.time = Sys.time()
number=500
temp3=dist_TF_cumsum_v2(data=All_Data_NAimputed[1:number,],id_annotation =unique(id_allAttributes[id_allAttributes$ids %in% 1:number,c("ids","Pwy")]) ,dist_metric = pcc_dist)
end.time = Sys.time()
end.time - start.time 

start.time = Sys.time()
number=500
temp3=dist_TF_cumsum_v3(data=All_Data_NAimputed[1:number,],attribute_list =attr_list(ECK.Pathway$ids,ECK.Pathway$Data) ,dist_metric = pcc_dist)
end.time = Sys.time()
end.time - start.time 

