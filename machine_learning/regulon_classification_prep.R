

#For each annotation set: how many members are there for a particualr annotation -> show by a barplot
id_regulon=id_allAttributes[,c("ids","regulator")] %>% unique
dim(id_regulon)

regulonAnnot=id_regulon$regulator[!is.na(id_regulon$regulator)]

regulonAnnotTable=sort(table(regulonAnnot),decreasing=T)

annotCounts=table(regulonAnnotTable)



#throw out annotations that appear only 1,2,3 times
annotCounts_filtered=annotCounts[as.numeric(names(annotCounts))>=4]

#barplot
barplot(annotCounts_filtered,xlab = "no. of Nichols' genes in the regulon'",ylab="no. of regulons")



##construct the label matrix of phenotypic profile to learn from


###filter out ids that are involved only in 1,2,3 gene pathways (or only 1,2,3 genes in Nichols are involved in those pathways)
regulonAnnotTable_atLeastFour=regulonAnnotTable[regulonAnnotTable>=4]
regulonToUse=names(regulonAnnotTable_atLeastFour)

id_regulon=id_allAttributes[,c("ids","regulator")] %>% unique()
attribute_list=attr_list(id_regulon$ids,id_regulon$regulator)


regulonLabel_df=as.data.frame(matrix(NA,ncol=length(regulonToUse),nrow=3979))

for(i in seq(regulonToUse)){
  
  regulon=regulonToUse[i]
  
  annotated_TFvec=sapply(attribute_list,FUN=function(annotations){
    regulon %in% annotations #If the "annotations" field is NA, this will return F
  })
  
  regulonLabel_df[,i]=annotated_TFvec
  
}

names(regulonLabel_df)=regulonToUse


save(regulonLabel_df,file="Data/regulonLabel_df.RData")
