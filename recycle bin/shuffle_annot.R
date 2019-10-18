#(!) The way the functions are written here might be the deadend (slow, too complicated to maintain)
#(!) Also I haven't finished the function shuffleAnnotation()


#Function that determines the total no. of difference between 2 "id_annot_list" obj. The ids and their order has to match.
diff_id_annot_list=function(id_annot_list,New){
  result=purrr::pmap(list(id_annot_list,New), function(annot,new_annot){
    identical(annot,as.character(NA))
    
    if(  !identical(annot,as.character(NA)) &  #If both of them are not NA
         !identical(new_annot,as.character(NA)) ){ 
      
      noOfDiff=sum(!(annot %in% new_annot)) + sum(!(new_annot %in% annot)) ##This is how I define no. of differences 
      
    }else if(identical(annot,as.character(NA)) &
             identical(new_annot,as.character(NA))
    ){ # If both of them are NA
      
      noOfDiff=sum(!(annot %in% new_annot)) + sum(!(new_annot %in% annot)) ##This is how I define no. of differences
      
    }else{ #If only 1 of them is NA
      
      noOfDiff=max(length(annot),length(new_annot))  
    }
    
    return(noOfDiff)
  }
  ) %>% unlist %>% sum
  
  return(result)
}

#Function that shuffles annotations in an "id_annot_list" obj
shuffleAnnotation=function(id_annot_list,fractionOfMisAnnot=c(0.005,0.05,0.25,0.5)){
  
  #total annotations (non-unique)
  tot=sapply(id_annot_list,FUN=function(annot){
    ifelse(!identical(annot,as.character(NA)),length(annot),0) #if there is no annotation (NA), the length should be 0
  }) %>% sum
  
  n=(tot*fractionOfMisAnnot) %>% round
  
  fake_annot_list=id_annot_list # will add all contaminated id_annot_list to it
  
  for(i in 1:length(n)){

      ##keep swaping annotations until the total difference matches the specified percentage
      New=id_annot_list
      count=1
      repeat{
        
        set.seed(i*count); ids=sample(names(New),2) #sample ids (first: where an annotation will be removed, second: where an annotation will be added)
        set.seed(i*count); annot=sample(id_annot_list[[ids[1]]],1) #sample annotation for that id
        
        
        New[[ids[1]]]=New[[ids[1]]][New[[ids[1]]]!=annot] #remove that annotation
        if(identical(New[[ids[1]]],character(0))) New[[ids[1]]]=NA #I would like the annotation slot to stay NA if it has been emptied
        
        New[[ids[2]]]=c(New[[ids[2]]],annot) #add that annotation to another id
        New[[ids[2]]]=New[[ids[2]]][!is.na(New[[ids[2]]])] #I have to manually remove NA if there were no annotation
        
        
        Diff=diff_id_annot_list(id_annot_list,New)
        print(Diff)  
        if(Diff>=n[i]){ 
          print(paste("Expected difference is",n[i]));  print(paste("real difference is",Diff)); break
        } 
        
        count=count+1
      }
      
    fake_annot_list[[i+1]]=New
  }

  
  return(fake_annot_list)
}


typeOfAnnotation="Pwy"
id_annot=id_allAttributes[,c("ids",typeOfAnnotation)] %>% unique
id_annot_list=attr_list(id_annot[[1]],id_annot[[2]])
str(id_annot_list)


test=shuffleAnnotation(id_annot_list,fractionOfMisAnnot=0.05)
