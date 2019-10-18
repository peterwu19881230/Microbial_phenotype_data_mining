#Make various annotation sets that have different levels of misannotations
##1.  A list that has original annotations + different no. of shuffled annotations
##2.  A list that has original annotations + different no. of annotations left out


class(id_allAttributes)
id_pwyAnnot=id_allAttributes[,c("ids","Pwy")] %>% unique
sum(!is.na(id_pwyAnnot$Pwy)) #2344 total pwy annotations (not unique)


##1.  A list that has original annotations + different no. of shuffled annotations

#The following mis-annotation step might not be the best (but the point is to somehow contaminate the database, and the following might be good enough)


#Remove ids that don't have any annotation (will add back after swaping annotations)
id_pwyAnnot_na_removed=id_pwyAnnot[!is.na(id_pwyAnnot$Pwy),]

fake_pwyAnnot=id_pwyAnnot_na_removed[,2:1]

#This creates no. of misannotations = 500, 1000, 1500, 2000
n=c(1,10,100,500,1000,1500,2000)
for(i in n){ #no. of id to be replaced (no. of misannotation to be done)
  
  set.seed(i); index_to_be_replaced=seq_along(fake_pwyAnnot$ids) %>% sample(i)
  
  contaminated_id=fake_pwyAnnot$ids
  for(index in index_to_be_replaced){
    set.seed(index) #Using index doesn't have any special meaning. Just to make the result reproducible
    
    pwy_for_this_index=fake_pwyAnnot$Pwy[index] #This should only have 1 (either it's a pwy or NA)
    
    if(is.na(pwy_for_this_index)){
      id_not_in_this_pwy=id_pwyAnnot$ids %>% unique
    }else{ 
      id_not_in_this_pwy=id_pwyAnnot$ids[id_pwyAnnot$Pwy!=pwy_for_this_index | is.na(id_pwyAnnot$Pwy)] %>% unique 
    } 
    
    
    #Change the id to a random id that is not the original id. 
    contaminated_id[index]=sample(id_not_in_this_pwy,1) 
    
  }
  
  fake_pwyAnnot[[paste("contaminated_ids_",i,sep="")]]=contaminated_id
  
  
  sum(id_pwyAnnot_na_removed$ids!=contaminated_id) %>% print #confirm the no. of misannotations
}



#add the ids that have no annotations back 
#(for different annotation sets in fake_pwyAnnot the ids to be added back are different)
fake_pwyAnnot_list=list()

for(i in 2:dim(fake_pwyAnnot)[2]){
  
  id_to_be_added=unique(id_pwyAnnot$ids)[!( unique(id_pwyAnnot$ids) %in% unique(fake_pwyAnnot[,i]) )]
  matrix_to_be_added=cbind(NA,id_to_be_added); colnames(matrix_to_be_added)=names(fake_pwyAnnot[,c(1,i)])
  
  fake_pwyAnnot_list[[i-1]]=rbind(fake_pwyAnnot[,c(1,i)],
                                  matrix_to_be_added)
} 


#Double check that they contain all the ids (1 ~ 3979)
map_dbl(fake_pwyAnnot_list,function(x){x[[2]] %>% unique %>% length}) 

#save(fake_pwyAnnot_list,file="Data/fake_pwyAnnot_list.RData")





##2.  A list that has original annotations + different no. of annotations left out
n=c(1,10,100,500,1000,1500,2000)

annotThrown_pwyAnnot=id_pwyAnnot
indexNotNA=(1:dim(annotThrown_pwyAnnot)[1])[!is.na(annotThrown_pwyAnnot$Pwy)]

for(i in seq(n)){
  
  set.seed(i)
  indexToBeThrown=sample(indexNotNA,n[i])
  
  new_annot=annotThrown_pwyAnnot$Pwy
  new_annot[indexToBeThrown]=NA
  
  annotThrown_pwyAnnot[[paste("less_pwy_",i,sep="")]]=new_annot
  
  #confirm the no. of annotations thrown out
  original=annotThrown_pwyAnnot$Pwy
  original[is.na(original)]="NA"
  
  later=new_annot
  later[is.na(later)]="NA"
  
  sum(original!=later) %>% print
}

#save(annotThrown_pwyAnnot,file="Data/annotThrown_pwyAnnot.RData")






