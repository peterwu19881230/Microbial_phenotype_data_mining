#(!)Be sure to consider No. of significant phenotypes for each condition
## convert all significant phenotypes to TRUE and others (including NAs) to FALSE 
significantData=(Binary_Data!=0 & !is.na(Binary_Data))
colSums(significantData) #No. of phenotypes for each condition


## Histogram of significant phenotypes for each condition. Minimum No. of phenotypes for a condition is 21
library(ggplot2)
qplot(main="Histogram of significant phenotypes for each condition",
      colSums(significantData),
      geom="histogram",
      bins = 15,
      col=I("black"),
      fill=I("white"),
      xlab="No. of phenotypes") 




# Test using threshold pcc=0.5
temp=throwPcc(sort.ConditionPCC.sql.NoIdent,0.5)

## How many unique conditions are left?
length(unique(c(as.character(temp[,1]),as.character(temp[,2])))) #140

##Create the sub-dataset that has only those 140 conditions
conditions=unique(c(as.character(temp[,1]),as.character(temp[,2])))
nonRedun.140cond.Data=All_Data[,conditions]
#save(nonRedun.140cond.Data,file="Data/sourced/nonRedun.140cond.Data.RData")
nonRedun.140cond.Ternary_Data=Ternary_Data[,conditions]
#save(nonRedun.140cond.Ternary_Data,file="Data/sourced/nonRedun.140cond.Ternary_Data.RData")

## pairwise pcc
nonRedun.140cond.cor_strain=cor(t(nonRedun.140cond.Data),use="pairwise.complete.obs",method="pearson")


#Reorder nonRedun.140cond.cor_strain into a SQL format
##Sort the output of a dist obj: https://stackoverflow.com/questions/31591546/sorting-the-output-of-dist
m1=nonRedun.140cond.cor_strain
m1[upper.tri(m1)]=0 
library(reshape2)
m2=subset(melt(m1), value!=0)
m2=m2[,c(2,1,3)] ##reorder the columns
names(m2)=c("strain1","strain2","Pearson.Correlation.Coefficient") ##name the columns
m2=m2[order(m2[,3],decreasing = T),] ##reorder by PCC 
sort.nonRedun.140cond.strainCC.sql.NoIdent=m2[m2$strain1!=m2$strain2,]
#save(sort.nonRedun.140cond.strainCC.sql.NoIdent,file="Data/sourced/sort.nonRedun.140cond.strainCC.sql.NoIdent.RData")



#Remodel the above so it looks like: 
#all cond -> throw 1 out, caculate the current highest pcc_cond -> store the index for the remaining cond + value of highest pcc -> repeat

throwPcc_nameSequence=throwPcc_sequence(sort.ConditionPCC.sql.NoIdent)

##This indices give the sequence of conditions that should be thrown out based on PCC (Each tim most correlated condition is thrown out)
indices=sapply(throwPcc_nameSequence,FUN = function(name){
  which( colnames(All_Data) %in% name)
})












