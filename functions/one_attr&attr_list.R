#Sub-functions required to run dist_TF_cumsum() and it's derivatives
#===============================================================================

# Functions to deal with relational data

## Compute the attribute list by inputting relational data
##name is the name vecotr, attribute is the vector of corresponding attributes
attr_list=function(name,attribute){
  uniqueName=unique(name)
  
  attribute_list=list()
  for(i in 1:length(uniqueName)){
    attribute_list[[i]]=attribute[name==uniqueName[i]]
  }
  
  names(attribute_list)=uniqueName
  
  return(attribute_list)
}

## Compute the vecotr of "at least 1 same attribute". NA is accepted for attribute.
one_attr=function(name,attribute){
  uniqueName=unique(name)
  
  combination=t(combn(uniqueName,2))
  
  attribute_list=list()
  for(i in 1:length(uniqueName)){
    TFvector=(name==uniqueName[i])
    attribute_list[[i]]=attribute[TFvector]
  }
  
  sameORnot=c()
  for(i in 1:dim(combination)[1]){
    id1=combination[i,1]
    id2=combination[i,2]
    
    sameORnot[i]=ifelse( sum( attribute_list[[id1]] %in% attribute_list[[id2]])>=1 & 
                           !(anyNA(attribute_list[[id1]])) & 
                           !(anyNA(attribute_list[[id2]])),
                         1,0)
  }
  
  return(data.frame(combination,sameORnot))
}