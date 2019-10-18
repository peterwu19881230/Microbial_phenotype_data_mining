shuffleAnnotation=function(id_annot,fractionOfMisAnnot=c(0.005,0.05,0.25,0.5)){
  tot=sum(!is.na(id_annot[[2]])) #no. of total annotations (not unique)
  
  #Remove ids that don't have any annotation (will add back after swaping annotations)
  id_annot_NAremoved=id_annot[!is.na(id_annot[[2]]),]
  
  fake_Annot=id_annot_NAremoved[,2:1] #I will swap the id to mis-annotate instead of swaping annotations
  
  #This creates no. of misannotations based on the total annotations
  n=(tot*fractionOfMisAnnot) %>% round
  
  for(i in n){ #no. of id to be replaced (no. of misannotation to be done)
    
    
    #set.seed(i) #I would add this back if I want a fixed result for a certain fraction of mis-annotations
    index_to_be_replaced=seq_along(fake_Annot[[2]]) %>% sample(i)
    
    contaminated_id=fake_Annot[[2]]
    for(index in index_to_be_replaced){
      set.seed(index) #Using index doesn't have any special meaning. Just to make the result reproducible
      
      annot_for_this_index=fake_Annot[[1]][index] #This should only have 1 (either it's a pwy or NA)
      
      if(is.na(annot_for_this_index)){
        id_not_in_this_annot=id_annot[[1]] %>% unique
      }else{ 
        id_not_in_this_annot=id_annot[[1]][id_annot[[2]]!=annot_for_this_index | is.na(id_annot[[2]])] %>% unique 
      } 
      
      
      #Change the id to a random id that is not the original id. 
      contaminated_id[index]=sample(id_not_in_this_annot,1) 
      
    }
    
    fake_Annot[[paste("contaminated_ids_",i,sep="")]]=contaminated_id
    
    
    sum(id_annot_NAremoved[[1]]!=contaminated_id) %>% print #confirm the no. of misannotations
  }
  
  
  #add the ids that have no annotations back 
  #(for different annotation sets in fake_Annot the ids to be added back are different)
  fake_annot_list=list()
  
  for(i in 2:dim(fake_Annot)[2]){
    
    id_to_be_added=unique(id_annot[[1]])[!( unique(id_annot[[1]]) %in% unique(fake_Annot[,i]) )]
    matrix_to_be_added=cbind(NA,id_to_be_added); colnames(matrix_to_be_added)=names(fake_Annot[,c(1,i)])
    
    resultDF=rbind(fake_Annot[,c(1,i)],
                   matrix_to_be_added)
    
    
    ##sort the dataframe to make the result cleaner
    resultDF=resultDF[order(resultDF[[2]] %>% as.numeric),]
    
    fake_annot_list[[i-1]]=resultDF
  }
  
  return(fake_annot_list)
}


#test using only "Pwy"
typeOfAnnotation="Pwy"
id_annot=id_allAttributes[,c("ids",typeOfAnnotation)] %>% unique


test=removeAnnotation(id_annot)

for(i in 1:5){
  identical(shuffleAnnotation(id_annot,fractionOfMisAnnot=rep(0.25,10)) ,
            shuffleAnnotation(id_annot,fractionOfMisAnnot=rep(0.25,10)) ) %>% print
}

