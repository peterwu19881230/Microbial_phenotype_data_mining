#New distance - cummulative sum function

##The "attribute_list" has to contain every id in the "data". Give "NA" if those ids don't have any attribute
##The data has to have corresponding ids from attribute_list
dist_TF_cumsum_matirxOperation=function(data,attribute_list,dist_metric){
  if(class(rownames(data))!="character")rownames(data)=as.character(rownames(data))
  
  #get table 1: id1_id2_sortedDistance
  #========================================
  distance=dist_metric(data)
  distance_table=melt_dist(distance)
  #========================================
  
  
  
  #get table 2: id1_id2_coAnnotated
  #========================================
  
  attribute_list=attribute_list[rownames(data)] 
  ##(important) change the order of elements in "attribute_list" so it matches those in "data"
  ##=>this prevents later problem in left_join. If this is not done the following won't be properly joined:
  ##=> Obj1 - Obj2 - dist = ("1","2","0.9"), Obj1 - Obj2 - sameORnot = ("2","1",1)
  
  
  uniqueAttr_vec=unlist(attribute_list) %>% unique
  uniqueAttr_vec=uniqueAttr_vec[!is.na(uniqueAttr_vec)]
  
  
  idRow_attrCol=sapply(attribute_list,FUN = function(attribute){ #idRow_attrCol will be created as a matrix
    uniqueAttr_vec %in% attribute  #c("A","B","C") %in% NA will give FALSE, so no worries here
      }) %>% t
  
  
  ##This gets pairwise comparison of having the same annotations or not
  coAnnotationMatrix=idRow_attrCol %*% t(idRow_attrCol)
  coAnnotationMatrix=ifelse(coAnnotationMatrix>=1,1,0) #The matrix structure will be preserved even after using ifelse()
  
  
  colnames(coAnnotationMatrix)=rownames(coAnnotationMatrix) #Have to manually give the colnames (which are identical to the row names)
  
  
  coAnnotated_table=coAnnotationMatrix %>% as.dist %>% melt_dist #after this, the result becomes a dataframe
  names(coAnnotated_table)[3]="sameORnot"  #have to correct the colname here
  #========================================
  
  #merge table 1 and table 2 and create a cumsum column
  
  result_df=cbind(distance_table,sameORnot=coAnnotated_table$sameORnot) 
  ##the order of the first 2 column should be the same. No need to use left_join
  ##(!) I am not sure if in rare cases this will not work 
  ##(eg. dist_metirc order the cols and rows in a different way -> melt_dist get different column 1 and column 2)
  
  result_df=result_df[order(result_df$Distance),]; rownames(result_df)=NULL #this resets the rownames after sorting
  result_df$cumsum=cumsum(result_df$sameORnot)
  
  return(result_df)
}





'

#Test on toy datasets:
#============================================================================================================
attribute_list=list(B=c("one", "two"),A=c("one", "two", "three", "four") , C=c("two", "four", "five", "six"), D=c("six", "seven"),E="four",f=NA,G=NA)
dat=data.frame(A=rnorm(1000),B=rnorm(1000),C=rnorm(1000),D=rnorm(1000),E=rnorm(1000),f=rnorm(1000),G=rnorm(1000)) %>% t


New=dist_TF_cumsum_matirxOperation(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
original=dist_TF_cumsum(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)

#I have to do some corrections because the old way:
##(i) returns a list 
##(ii) the name of the distance column is "distance" Instead of "Distance"
##(iii) the rownames are messed up
original=original[[1]]; names(original)[3]="Distance"; rownames(original)=as.character(1:dim(original)[1])

#identical(New,original) #might be some nuances that cause this to be FALSE
#all.equal(New,original)

#This confirms that all columns are the same
mapply(New,original,FUN = identical) #mapply goes columnwise (because dataframe is a special case of list with columns as elements) 
#============================================================================================================



#Test on real datasets:
#============================================================================================================
dat=All_Data_NAimputed
id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique()
attribute_list=attr_list(id_pwy$ids,id_pwy$Pwy)


start.time = Sys.time()
original=dist_TF_cumsum(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
end.time = Sys.time()
end.time - start.time #Time difference of 9.092783 mins


start.time = Sys.time()
New=dist_TF_cumsum_matirxOperation(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
end.time = Sys.time()
end.time - start.time #Time difference of 39.23815 secs. More than 10 times faster than the original method


#I have to do some corrections because the old way:
##(i) returns a list 
##(ii) the name of the distance column is "distance" Instead of "Distance"
##(iii) the rownames are messed up
original=original[[1]]; names(original)[3]="Distance"; rownames(original)=as.character(1:dim(original)[1])

#identical(original,New) #might be some nuances that cause this to be FALSE
#all.equal(New,original)

#This confirms that all columns are the same
mapply(New,original,FUN = identical) #mapply goes columnwise (because dataframe is a special case of list with columns as elements) 
#============================================================================================================

'




