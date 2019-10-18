

#For each annotation set: how many members are there for a particualr annotation -> show by a barplot
id_operon=id_allAttributes[,c("ids","operon")] %>% unique
dim(id_operon)

operonAnnot=id_operon$operon[!is.na(id_operon$operon)]

operonAnnotTable=sort(table(operonAnnot),decreasing=T)

annotCounts=table(operonAnnotTable)



#throw out annotations that appear only 1,2,3 times
annotCounts_filtered=annotCounts[as.numeric(names(annotCounts))>=4]

#barplot
barplot(annotCounts_filtered,xlab = "no. of Nichols' genes in the operon'",ylab="no. of operons")



##construct the label matrix of phenotypic profile to learn from


###filter out ids that are involved only in 1,2,3 gene operons (or only 1,2,3 genes in Nichols are involved in those operons)
operonAnnotTable_atLeastFour=operonAnnotTable[operonAnnotTable>=4]
operonToUse=names(operonAnnotTable_atLeastFour)

id_operon=id_allAttributes[,c("ids","operon")] %>% unique()
attribute_list=attr_list(id_operon$ids,id_operon$operon)


operonLabel_df=as.data.frame(matrix(NA,ncol=length(operonToUse),nrow=3979))

for(i in seq(operonToUse)){
  
  operon=operonToUse[i]
  
  annotated_TFvec=sapply(attribute_list,FUN=function(annotations){
    operon %in% annotations #If the "annotations" field is NA, this will return F
  })
  
  operonLabel_df[,i]=annotated_TFvec
  
}

names(operonLabel_df)=operonToUse


save(operonLabel_df,file="Data/operonLabel_df.RData")