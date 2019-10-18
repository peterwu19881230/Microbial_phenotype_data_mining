#Distribution of GO annotations made to Nichols'
##1. using all the strains (BP, MF, CC are separated)
##2. using pair of strains that have similaity > 95% of the pairs (BP, MF, CC are separated)

id_GO=id_allAttributes[,c("ids","GO")] %>% unique
id_GO_noNA=id_GO[!is.na(id_GO$GO),]
dim(id_GO_noNA)


id_GO_BP=id_GO_noNA$GO[AnnotationDbi::Ontology(id_GO_noNA$GO)=="BP" & !is.na(AnnotationDbi::Ontology(id_GO_noNA$GO))]
id_GO_MF=id_GO_noNA$GO[AnnotationDbi::Ontology(id_GO_noNA$GO)=="MF" & !is.na(AnnotationDbi::Ontology(id_GO_noNA$GO))]
id_GO_CC=id_GO_noNA$GO[AnnotationDbi::Ontology(id_GO_noNA$GO)=="CC" & !is.na(AnnotationDbi::Ontology(id_GO_noNA$GO))]


GO_BP_freq=table(id_GO_BP) %>% sort(decreasing=T)
barplot(GO_BP_freq[1:100],las = 2,cex.names=0.8) #see top 100
AnnotationDbi::Term(names(GO_BP_freq[1:10])) #description for the top 10 terms
length(unique(id_GO_BP)) #total no. of terms

GO_MF_freq=table(id_GO_MF) %>% sort(decreasing=T)
barplot(GO_MF_freq[1:100],las = 2,cex.names=0.8) #see top 100
AnnotationDbi::Term(names(GO_MF_freq[1:10])) #description for the top 10 terms
length(unique(id_GO_MF)) #total no. of terms

GO_CC_freq=table(id_GO_CC) %>% sort(decreasing=T)
barplot(GO_CC_freq,las = 2,cex.names=0.8) #total no.=139. I just plot everything
AnnotationDbi::Term(names(GO_CC_freq[1:10])) #description for the top 10 terms
length(unique(id_GO_CC)) #total no. of terms




#What terms are for the strains that have similarity > 95% of all the strain pairs (haven't finished)


##BP==================================================================================================
orderedBy_BP=strain1_strain2_pcc_WangBP_noNA[order(strain1_strain2_pcc_WangBP_noNA$BP,decreasing=T),]
dim(orderedBy_BP) #3465028*0.05 ~ 173251 (similarity=0.547)

#Get the strain pairs that have > 0.547 similarity
BPtop5_strains=c(orderedBy_BP[1:173251,"strain1"],orderedBy_BP[1:173251,"strain2"]) #I want to do frequency count for GO terms so I don't use unique()


#Filter GO annotations so there are only BP terms (be careful about obsoleted terms)
id_GO_noNA_BP=id_GO_noNA[AnnotationDbi::Ontology(id_GO_noNA$GO)=="BP" &
                           !is.na(AnnotationDbi::Ontology(id_GO_noNA$GO)),] 
##Some terms are obsolete (eg. GO:0022891, GO:0046909) and AnnotationDbi::Ontology() will return NA, so I have to deal with that situation


#Get all terms for all the above strains
BP_terms_forStrains=sapply(as.character(BPtop5_strains),FUN = function(strain){
  id_GO_noNA_BP$GO[id_GO_noNA_BP$ids==strain]
})

str(BP_terms_forStrains)

BP_terms_forStrains_unlisted=unlist(BP_terms_forStrains)

BP_dat=table(BP_terms_forStrains_unlisted) %>% sort(decreasing=T)

length(BP_dat) #1395

#barplot of the first 100
barplot(BP_dat[1:100],las = 2,cex.names=0.8)

#Get description for the top 10 terms
AnnotationDbi::Term(names(BP_dat[1:10]))
##==================================================================================================





##MF==================================================================================================
orderedBy_MF=strain1_strain2_pcc_WangMF_noNA[order(strain1_strain2_pcc_WangMF_noNA$MF,decreasing=T),]
dim(orderedBy_MF) #3086370*0.05 ~ 154319 (similarity=0.815)

#Get the strain pairs that have > 0.547 similarity
MFtop5_strains=c(orderedBy_MF[1:154319,"strain1"],orderedBy_MF[1:154319,"strain2"]) #I want to do frequency count for GO terms so I don't use unique()


#Filter GO annotations so there are only MF terms (obso)
id_GO_noNA_MF=id_GO_noNA[AnnotationDbi::Ontology(id_GO_noNA$GO)=="MF" &
                           !is.na(AnnotationDbi::Ontology(id_GO_noNA$GO)),] 
##It's possible that some terms are obsolete and AnnotationDbi::Ontology() will return NA, so I have to deal with that situation


#Get all terms for all the above strains
MF_terms_forStrains=sapply(as.character(MFtop5_strains),FUN = function(strain){
  id_GO_noNA_MF$GO[id_GO_noNA_MF$ids==strain]
})

str(MF_terms_forStrains)

MF_terms_forStrains_unlisted=unlist(MF_terms_forStrains)

MF_dat=table(MF_terms_forStrains_unlisted) %>% sort(decreasing=T)

length(MF_dat) #1474

#barplot of the first 100
barplot(MF_dat[1:100],las = 2,cex.names=0.8)

#Get description for the top 10 terms
AnnotationDbi::Term(names(MF_dat[1:10]))
##==================================================================================================





#CC==================================================================================================
orderedBy_CC=strain1_strain2_pcc_WangCC_noNA[order(strain1_strain2_pcc_WangCC_noNA$CC,decreasing=T),]
dim(orderedBy_CC) #2579856*0.05 ~ 128992 (similarity=1)

#Get the strain pairs that have > 0.547 similarity
CCtop5_strains=c(orderedBy_CC[1:154442,"strain1"],orderedBy_CC[1:154442,"strain2"]) #I want to do frequency count for GO terms so I don't use unique()


#Filter GO annotations so there are only CC terms (obso)
id_GO_noNA_CC=id_GO_noNA[AnnotationDbi::Ontology(id_GO_noNA$GO)=="CC" &
                           !is.na(AnnotationDbi::Ontology(id_GO_noNA$GO)),] 
##It's possible that some terms are obsolete and AnnotationDbi::Ontology() will return NA, so I have to deal with that situation


#Get all terms for all the above strains
CC_terms_forStrains=sapply(as.character(CCtop5_strains),FUN = function(strain){
  id_GO_noNA_CC$GO[id_GO_noNA_CC$ids==strain]
})

str(CC_terms_forStrains)

CC_terms_forStrains_unlisted=unlist(CC_terms_forStrains)

CC_dat=table(CC_terms_forStrains_unlisted) %>% sort(decreasing=T)

length(CC_dat) 

#barplot of all (total no.=124)
barplot(CC_dat,las = 2,cex.names=0.8)

#Get description for the top 10 terms
AnnotationDbi::Term(names(CC_dat[1:10]))
##==================================================================================================