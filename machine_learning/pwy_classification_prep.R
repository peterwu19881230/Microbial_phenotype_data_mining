

#For each annotation set: how many members are there for a particualr annotation -> show by a barplot
id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique
dim(id_pwy)

pwyAnnot=id_pwy$Pwy[!is.na(id_pwy$Pwy)]

pwyAnnotTable=sort(table(pwyAnnot),decreasing=T)


annotCounts=c()
for(i in unique(pwyAnnotTable)){
  annotCounts=c(annotCounts,sum(pwyAnnotTable==i))
}
names(annotCounts)=unique(pwyAnnotTable)

print(annotCounts)



#throw out annotations that appear only 1,2,3 times
annotCounts_filtered=annotCounts[as.numeric(names(annotCounts))>=4]

#barplot
barplot(annotCounts_filtered,xlab = "no. of Nichols' genes in the pathway'",ylab="no. of pathways")


#for genes that have more than 1 pathway, maybe I should 
#1. randomly delete annotations until there's only 1 annotation to that gene -> but this way genes involved in multiple pathways that have similar profile won't be catched
#2. propogate the phenotypic profiles (let some profiles appear more than once in the dataset to be trained) -> This will result in having identical profiles that have different labels. Would this be chaotic?
#3. create a label matrix, where there are true-false labels for each pathway


#try 2. :


##construct the phenotypic profile to learn from

###filter out ids that are involved only in 1,2,3 gene pathways (or only 1,2,3 genes in Nichols are involved in those pathways)
pwyAnnotTable_atLeastFour=pwyAnnotTable[pwyAnnotTable>=4]
pwyToUse=names(pwyAnnotTable_atLeastFour)

idToUse=id_pwy$ids[id_pwy$Pwy %in% pwyToUse] #This is a character vector, but rownames(All_Data_NAimputed) is also character

#How many duplicated ids (resulting in duplicated profiles) are there:
#duplicatedID=table(table(idToUse))
#barplot(duplicatedID)


pwyNewPheno=All_Data_NAimputed[idToUse,]
dim(pwyNewPheno)[1] #1954 

pwyLabel=id_pwy$Pwy[id_pwy$Pwy %in% pwyToUse]

save(pwyNewPheno,pwyLabel,file="Data/phenoForPythonMachineLearning.RData")

#It turned out that duplicating the profiles resulted in chaotic result. See: pwy_classification.py



#try 3. :

##construct the label matrix of phenotypic profile to learn from


###filter out ids that are involved only in 1,2,3 gene pathways (or only 1,2,3 genes in Nichols are involved in those pathways)
pwyAnnotTable_atLeastFour=pwyAnnotTable[pwyAnnotTable>=4]
pwyToUse=names(pwyAnnotTable_atLeastFour)

id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique()
attribute_list=attr_list(id_pwy$ids,id_pwy$Pwy)


pwyLabel_df=as.data.frame(matrix(NA,ncol=length(pwyToUse),nrow=3979))

for(i in seq(pwyToUse)){
  
  pwy=pwyToUse[i]
  
  annotated_TFvec=sapply(attribute_list,FUN=function(annotations){
    pwy %in% annotations #If the "annotations" field is NA, this will return F
  })
  
  pwyLabel_df[,i]=annotated_TFvec
  
}
  
names(pwyLabel_df)=pwyToUse


save(pwyLabel_df,file="Data/pwyLabel_df.RData")








