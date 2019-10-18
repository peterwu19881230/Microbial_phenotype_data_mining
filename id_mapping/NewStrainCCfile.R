#Make a strainCC file with all the gene names synonyms (lots of them will have different names but the same data)

#A new strainCC object
strainCC=cor(t(All_Data[,-1]),use="pairwise.complete.obs",method="pearson") ##This is the question for Nasos (whether they used pairwise.complete.obs)



##rename the genes that have NULL (total=80) for EcoCyc ID to be NULL1, NULL2, NULL3...
new.names=c()
j=0
for(i in 1:length(ECK_1st_table$EcoCycID)){
  if(ECK_1st_table$EcoCycID[i]=="NULL"){
    j=j+1
    new.names[i]=paste("NULL",j,sep="")
  }else{new.names[i]=ECK_1st_table$EcoCycID[i]}
}

temp=cbind(new.names,strainCC)
temp=rbind(c("/",new.names),new.strainCC)
write.table(temp,"new.strainCC.txt",row.names = F,col.names = F,quote = F)
rm(temp)





