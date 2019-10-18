#Rcpp practice: Using bottle neck of my self-written function 


##Synthetic data
library(dplyr)
distance_table=t(combn(c("A","B","C","D","E"),2)) %>% as.data.frame
attribute_list=list(A=c("a","b","c"),B=NA,C=c("a","b"),D=c("e","f"),E=c("b","e","g"))


#original R code
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
#====================================================================================



#Rcpp
library(Rcpp)
#setwd("/Users/peterwu/Google Drive File Stream/Team Drives/HuLab/Nichols_Data_mining/")
sourceCpp("test/Rcpp_practice.cpp")
sourceCpp("test/coannoted_or_not.cpp")


# Compare speeds of 2 ways:

##run the toy example many times:
n=100000

start.time = Sys.time()
for(i in 1:n){
  #original R code
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
  
  #==================================================================================== 
}
end.time = Sys.time()
end.time - start.time #Time difference of 1.312657 mins


start.time = Sys.time()
for(i in 1:n){
  coannoted_or_not(distance_table,attribute_list) 
}
end.time = Sys.time()
end.time - start.time #Time difference of 12.43087 secs



##run larger samples:

distance_table=strain1_strain2
distance_table$strain1=as.character(distance_table$strain1)
distance_table$strain2=as.character(distance_table$strain2)

id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique()
attribute_list=attr_list(id_pwy$ids,id_pwy$Pwy)

start.time = Sys.time()
  #original R code
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
  
  #==================================================================================== 
end.time = Sys.time()
end.time - start.time #Time difference of 7.120816 mins




start.time = Sys.time()

result=coannoted_or_not(distance_table,attribute_list) 

end.time = Sys.time()
end.time - start.time  #Time difference of 7.23962 mins


#Rcpp version produces the same result
identical(result,as.integer(sameORnot))



start.time = Sys.time()

result2=coannoted_or_not_parallel(distance_table,attribute_list,7)  #use 7 cores

end.time = Sys.time()
end.time - start.time  #Time difference of 6.950679 mins on PC. Not faster. Weird.





