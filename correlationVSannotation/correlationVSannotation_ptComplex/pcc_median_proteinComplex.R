#Goal: calculate the average (or median) pcc of genotype pcc pairs that are in the same protein complex

#Get the pcc matrix that contains only strains involved in protein complexes
##Use ids.originalName.pcomplex defined in proteinComplex.R

pt.list=split(ids.originalName.pcomplex,ids.originalName.pcomplex$pcomplex)

ptcomplex.pcc.matrix=sapply(pt.list,FUN=function(x){
  strains=x$id
  index.matrix=t(combn(strains,2)) ##combn() creates a matrix that contains all pairwise combination of strains
  
  pccs=c()
  for(i in 1:dim(index.matrix)[1]){
    pccs=c(pccs,cor_strains[index.matrix[i,1],index.matrix[i,2]]) 
    ##The order of index.matrix[i,1] and index.matrix[i,2] doesn't matter. Should get the same pcc
  }
  avg.pcc=mean(pccs)
  median.pcc=median(pccs)
  return(c(avg.pcc=avg.pcc,median.pcc=median.pcc,NoOfGenesUsed=dim(x)[1]))
}
)

#Matrix -> (transpose) -> Dataframe conversion
ptcomplex.pcc.dataframe=as.data.frame(t(ptcomplex.pcc.matrix))


#Distribution of mean pccs
hist(ptcomplex.pcc.dataframe$avg.pcc,main="Distribution of mean pccs of protein complexes",xlab="Pearson Correlation Coefficient")
summary(ptcomplex.pcc.dataframe$avg.pcc)

#Distribution of median pccs
hist(ptcomplex.pcc.dataframe$median.pcc,main="Distribution of median pccs of protein complexes",xlab="Pearson Correlation Coefficient")
summary(ptcomplex.pcc.dataframe$median.pcc)

#Negative control (average of all pcc)
mean(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient) ##0.0007193512



