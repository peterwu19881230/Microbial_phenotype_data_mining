

#For each annotation set: how many members are there for a particualr annotation -> show by a barplot
id_pcomplex=id_allAttributes[,c("ids","pcomplex")] %>% unique
dim(id_pcomplex)

pcomplexAnnot=id_pcomplex$pcomplex[!is.na(id_pcomplex$pcomplex)]

pcomplexAnnotTable=sort(table(pcomplexAnnot),decreasing=T)

annotCounts=table(pcomplexAnnotTable)



#throw out annotations that appear only 1,2,3 times
annotCounts_filtered=annotCounts[as.numeric(names(annotCounts))>=4]

#barplot
barplot(annotCounts_filtered,xlab = "no. of Nichols' genes in the pcomplex'",ylab="no. of pcomplexs")



##construct the label matrix of phenotypic profile to learn from


###filter out ids that are involved only in 1,2,3 gene pathways (or only 1,2,3 genes in Nichols are involved in those pathways)
pcomplexAnnotTable_atLeastFour=pcomplexAnnotTable[pcomplexAnnotTable>=4]
pcomplexToUse=names(pcomplexAnnotTable_atLeastFour)

id_pcomplex=id_allAttributes[,c("ids","pcomplex")] %>% unique()
attribute_list=attr_list(id_pcomplex$ids,id_pcomplex$pcomplex)


pcomplexLabel_df=as.data.frame(matrix(NA,ncol=length(pcomplexToUse),nrow=3979))

for(i in seq(pcomplexToUse)){
  
  pcomplex=pcomplexToUse[i]
  
  annotated_TFvec=sapply(attribute_list,FUN=function(annotations){
    pcomplex %in% annotations #If the "annotations" field is NA, this will return F
  })
  
  pcomplexLabel_df[,i]=annotated_TFvec
  
}

names(pcomplexLabel_df)=pcomplexToUse


save(pcomplexLabel_df,file="Data/pcomplexLabel_df.RData")
