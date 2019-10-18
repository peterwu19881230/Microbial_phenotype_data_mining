## Compute the vecotr of "at least 1 same attribute". NA is accepted for attribute.


# This is a much faster implementation
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


#This is the same function written much neater. But it's terribly slow
one_attr=function(name,attribute){
  combination=t(combn(unique(name),2))
  
  sameORnot=c()
  for(i in 1:dim(combination)[1]){
    name1=combination[i,1]
    name2=combination[i,2]
    
    sameORnot[i]=ifelse( sum( attribute[name==name1] %in% attribute[name==name2])>=1 & 
                           !(anyNA(attribute[name==name1])) & 
                           !(anyNA(attribute[name==name2])),
                         1,0)
  }
  
  return(data.frame(combination,sameORnot))
}


## Test using real data
dat=ECK.Pathway
dat$Data[dat$Data=="Not in any pathway"]=NA

start=Sys.time()
vec=one_attr(dat$ids[1:1000],dat$Data[1:1000])
end=Sys.time()
end-start 
