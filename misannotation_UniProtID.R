#generate a dataframe where UniProtID is shuffled: 
# => NicholsId - EcoCycID - shuffled UniProtID 

id_EcoCyc_UniProt=id_allAttributes[,c("ids","EcoCycID","UniProtID")] %>% unique
dim(id_EcoCyc_UniProt)

#Check the distribution of UniProtIDs
##id_EcoCyc_UniProt$UniProtID %>% table %>% sort(decreasing=T)


id_EcoCyc_originalUniProt_misUniProt=id_EcoCyc_UniProt

set.seed(102)
id_EcoCyc_originalUniProt_misUniProt$randomized=sample(id_EcoCyc_UniProt$UniProtID,
                                                       length(id_EcoCyc_UniProt$UniProtID))


write.table(id_EcoCyc_originalUniProt_misUniProt,file="Data/id_EcoCyc_originalUniProt_misUniProt.txt",
            row.names=F,col.names=T,quote=F,sep = "\t") 
