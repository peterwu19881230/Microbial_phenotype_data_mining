#Can know phenotype signature determine whether a gene (or genotype) is involved in a certain pathway?





# one of the element in PWY-901: Nichols' id=210
test_id=1481
phenoSig=as.numeric(All_Data[test_id,])

pcc.phenoSig_to_pwyPCC=list()
for(i in 1:length(unique_pwy_annot)){
  ids=pathway_ids[[unique_pwy_annot[i]]]
  corr=c()
  for(j in 1:length(ids)){
    if(ids[j]==test_id){next} #Skip when it's correlating to itself
    aGeneInThatPwy=as.numeric(All_Data[ids[j],])
    corr=c(corr,cor(phenoSig,aGeneInThatPwy,use="pairwise.complete.obs",method="pearson")) # can modify this: can just call the numbers from cor_strain
  }
  pcc.phenoSig_to_pwyPCC[[i]]=mean(corr) #Can try median() later
}


#Convert to dataframe to simplify the output
pcc.phenoSig_to_pwyPCC=cbind(unique_pwy_annot,as.numeric(pcc.phenoSig_to_pwyPCC))





#Some stuff that might be useful:
#in_pwy_table$ids=as.numeric(in_pwy_table$ids)
#unique_id=unique(in_pwy_table$ids)
#Pwy_Data=All_Data[unique_id,]
#in_pwy_cor_strains=cor(t(Pwy_Data),use="pairwise.complete.obs",method="pearson") ##This is the question for Nasos (whether they used pairwise.complete.obs)

