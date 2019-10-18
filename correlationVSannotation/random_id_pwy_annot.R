#The following is my own way of defining completely random annotation. Not sure if it's correct

pwy_annot=id_pwyAnnot$Pwy[!is.na(id_pwyAnnot$Pwy)] 
##I don't take the unique annotations since I want the no. of each 
##unique annotations to stay the same

uniq_id=unique(id_pwyAnnot$ids) #There is no NA in id, so no need to remove it

## 1. For each pwy_annot assign a random id 
## 2. The rest of the ids that don't have annotations get an NA annotation
new_id=c()
for(i in seq_along(pwy_annot)){
  set.seed(i)
  new_id[i]=sample(uniq_id,1)
}

id_not_assigned=uniq_id[!(uniq_id %in% new_id)]
length(unique(new_id))+length(id_not_assigned) # Double check. This should give 3979

bad_id_pwy_annotations=cbind(c(new_id,id_not_assigned),
      c(pwy_annot,rep(NA,length(id_not_assigned))))


##This is an older way of defining random that I have thought about
##bad_pwy_annotations=transmute(fake_pwyAnnot,ids=sample(ids))


