##1. Function to generate a list that has original annotations + different no. of shuffled annotations
##2. Function to generate a list that has original annotations + different no. of annotations left out
##3. Function to create mis-annotation set based on types of annotations
##4. #Function to create strain1_strain2_coAnnotation
#(modified from contaminate_annotation.R)

##id_annot: df that has all pairs of id-annotations (including ids that have no-annotations)
##fractionOfMisAnnot is a numeric vector that decides the ratio of no. of mis-annotations being left out
##(Question)Should I add random? How do I even define random?


##1. Function to generate a list that has original annotations + different no. of shuffled annotations
#================================================================================================
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
    
    ids=fake_Annot[[2]]
    for(index in index_to_be_replaced){
      set.seed(index) #Using index doesn't have any special meaning. Just to make the result reproducible
      
      annot_for_this_index=fake_Annot[[1]][index] #This should only have 1 (either it's a pwy or NA)
      
      if(is.na(annot_for_this_index)){
        id_not_in_this_annot=id_annot[[1]] %>% unique
      }else{ 
        id_not_in_this_annot=id_annot[[1]][id_annot[[2]]!=annot_for_this_index | is.na(id_annot[[2]])] %>% unique 
      } 
      
      
      #Change the id to a random id that is not the original id. 
      ids[index]=sample(id_not_in_this_annot,1) 
      
    }
    
    fake_Annot=cbind(fake_Annot,ids)
    
    #Can add this line back if needed:
    #sum(id_annot_NAremoved[[1]]!=ids) %>% print #confirm the no. of misannotations
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
#================================================================================================

#2. Function to generate a list that has original annotations + different no. of annotations left out
#================================================================================================
removeAnnotation=function(id_annot,fractionOfMisAnnot=c(0.005,0.05,0.25,0.5)){
  
  tot=sum(!is.na(id_annot[[2]])) #no. of total annotations (not unique)
  
  
  #This creates no. of misannotations based on the total annotations
  n=(tot*fractionOfMisAnnot) %>% round
  
  
  annot_id=id_annot[,2:1]
  annot_id=annot_id[order(annot_id[[2]] %>% as.numeric),] ##sort the original dataframe to make the result cleaner
  
  fake_annot_list=list(annot_id) #attach the original id_annot and then add the contaminated ones
  count=2
  for(i in n){ 
    
    nonNA_index=which(!is.na(annot_id[[1]]))
    
    # set.seed(i) #I would add this back if I want a fixed result for a certain fraction of mis-annotations
    index=sample(nonNA_index,i) #sample an annotation slot that is not NA
    
    annot_new_id=annot_id
    annot_new_id[[1]][index]=NA #remove the annotations (set them to NA)
    
    
    #Can add this line back if needed:
    #Diff=dim(annot_id)[1] - sum( mapply(identical,annot_id[[1]],annot_new_id[[1]]) );  print(Diff) # print #confirm the no. of misannotations
    
    
    #Remove redundant ids (ids that have at least 1 annotation while having additional NA spots. I will remove those NA spots)
    annot_new_id=annot_new_id[!is.na(annot_new_id[[1]]),] #remove all indices that are NA
    id_to_bind=unique(annot_id[[2]])[!( unique(annot_id[[2]]) %in% annot_new_id[[2]] )]
    annot_id_to_bind=cbind(NA,id_to_bind); colnames(annot_id_to_bind)=colnames(annot_new_id) #(Note: annot_id_to_bind is a matrix)
    annot_new_id=rbind(annot_new_id,annot_id_to_bind) #bind id with no annotation back #(Note: annot_new_id is a dataframe)
    
    
    ##sort the dataframe to make the result cleaner
    annot_new_id=annot_new_id[order(annot_new_id[[2]] %>% as.numeric),]
    
    
    fake_annot_list[[count]]=annot_new_id
    count=count+1
  }
  
  return(fake_annot_list)
}
#================================================================================================


##test example
#class(id_allAttributes)
#id_pwyAnnot=id_allAttributes[,c("ids","Pwy")] %>% unique



#fake_pwyAnnot_list=shuffleAnnotation(id_pwyAnnot)
#map_dbl(fake_pwyAnnot_list,function(x){x[[2]] %>% unique %>% length}) ##Double check that they contain all the ids (1 ~ 3979)

#fake_pwyAnnot_list=removeAnnotation(id_pwyAnnot)
#map_dbl(fake_pwyAnnot_list,function(x){x[[2]] %>% unique %>% length}) ##Double check that they contain all the ids (1 ~ 3979)


#Function to create mis-annotation set based on types of annotations
#================================================================================================
createMisAnnotSet=function(typeOfAnnotation){
  id_Annot=id_allAttributes[,c("ids",typeOfAnnotation)] %>% unique
  
  shuffledAnnot_list=shuffleAnnotation(id_Annot)
  removedAnnot_list=removeAnnotation(id_Annot)
  
  return(list(shuffledAnnot_list=shuffledAnnot_list,
              removedAnnot_list=removedAnnot_list))
}
#================================================================================================



#Function to create strain1_strain2_coAnnotation
#================================================================================================
pairwiseCoannotation=function(annot_id){
  
  attribute_list=attr_list(annot_id[[2]],annot_id[[1]])
  
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
  
  return(coAnnotated_table)
}
#================================================================================================





