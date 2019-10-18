#Change the EcoGene ID at the bottom of this file to UniProt ID: ecocyc.gaf
##This houses the GAF files:
##http://current.geneontology.org/products/pages/downloads.html
## Ours (E. coli) will be listed under multi-species and will be ecocyc.gaf


ecocyc.gaf=read.table("Data/ecocyc.2.12.2019.gaf",sep="\t",quote="",comment.char="!") #quote="" is to disable quoting

names(id_allAttributes)

EcoCycID_UniProtID=unique(id_allAttributes[,c("EcoCycID","UniProtID")])
clean_EcoCycID_UniProtID=EcoCycID_UniProtID[EcoCycID_UniProtID[,1]!="NULL" &
                                            !is.na(EcoCycID_UniProtID[,2])  ,] #remove names that don't have mapping
dim(clean_EcoCycID_UniProtID)

hasUniProt=(ecocyc.gaf[[1]]=="UniProtKB")
ecocyc.gaf_keep=ecocyc.gaf[hasUniProt,]
ecocyc.gaf_to_be_switched=ecocyc.gaf[!hasUniProt,]


#Double check that there is no duplications for either column
clean_EcoCycID_UniProtID[,"EcoCycID"][clean_EcoCycID_UniProtID[,"EcoCycID"] %>% duplicated]
clean_EcoCycID_UniProtID[,"UniProtID"][clean_EcoCycID_UniProtID[,"UniProtID"] %>% duplicated]
##=> No duplications. One-to-One relationship


#Finish the following after the above problem is solved:

#The 2nd column is to be switched (EcoGene ID -> UniProt ID)

mapped_uniprotID=sapply(ecocyc.gaf_to_be_switched[[2]],function(EcoGene){
  uniprotID=clean_EcoCycID_UniProtID$UniProtID[clean_EcoCycID_UniProtID$EcoCycID==EcoGene]
  
  if(identical(uniprotID,character(0))){
    return(NA)
  }else return(uniprotID) 
})


ecocyc.gaf_hasBeenSwitched=ecocyc.gaf_to_be_switched
ecocyc.gaf_hasBeenSwitched[[1]]=rep("UniProtKB",length(ecocyc.gaf_hasBeenSwitched[[1]]))
ecocyc.gaf_hasBeenSwitched[[2]]=mapped_uniprotID


ecocyc.gaf_EcoGene_switchedTo_UniProt=rbind(ecocyc.gaf_keep,ecocyc.gaf_hasBeenSwitched)

#double check that the dimension of the new obj is right
dim(ecocyc.gaf_EcoGene_switchedTo_UniProt)
dim(ecocyc.gaf)

write.table(ecocyc.gaf_EcoGene_switchedTo_UniProt,file="Data/ecocyc_EcoGene_switchedTo_UniProt.gaf",row.names=F,col.names = F,quote = F, sep = "\t")



