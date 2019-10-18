#Pairwise semantic similarity of Nichols strains (by sets of GO annotations (BP, MF, CC are applied separately) )




id_allGO=id_allAttributes[,c("ids","GO")] %>% unique
id_allGO=id_allGO[!is.na(id_allGO$GO),]

id_BP=id_allGO[AnnotationDbi::Ontology(id_allGO$GO)=="BP" & !is.na(AnnotationDbi::Ontology(id_allGO$GO)),]
id_MF=id_allGO[AnnotationDbi::Ontology(id_allGO$GO)=="MF" & !is.na(AnnotationDbi::Ontology(id_allGO$GO)),]
id_CC=id_allGO[AnnotationDbi::Ontology(id_allGO$GO)=="CC" & !is.na(AnnotationDbi::Ontology(id_allGO$GO)),]



pairwise_GOannotation=function(id_GO){
  strain1strain2_transposed=as.data.frame((combn(1:3979,2))) 
  
  library(parallel)
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  
  
  clusterExport(cl=cl,"id_GO",envir=environment()) 
  #envir=environment() allows clusterExport to search the environment within this function
  #ref: https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function/12024448
  
  
  GOannotation=parLapply(cl=cl,strain1strain2_transposed,fun=function(strain1strain2){ 
    ##Although input is not a list, it will be coerced by as.list() according to ?lapply
    
    strain1=as.character(strain1strain2[1])
    strain2=as.character(strain1strain2[2])
    
    annotation1=id_GO$GO[id_GO$ids==strain1]
    annotation2=id_GO$GO[id_GO$ids==strain2]
    
    return(list(annotation1,annotation2))
    
  })
  
  
  stopCluster(cl)
  return(GOannotation)
}


if(!file.exists("Data/sourced/Wang_pairwise_similarity.RData")){

start.time=Sys.time()

pairwise_id_BP=pairwise_GOannotation(id_BP)
pairwise_id_MF=pairwise_GOannotation(id_MF)
pairwise_id_CC= pairwise_GOannotation(id_CC)

#Note: id without GO annotations have character(0) in the list
save(pairwise_id_BP,pairwise_id_MF,pairwise_id_CC,file="Data/sourced/pairwise_GOAnnotation_separated.RData")


end.time=Sys.time()
end.time-start.time #Time difference of 16.3797 mins





##parallelized code to get distances from Wang method 
#BiocManager::install("AnnotationHub", version = "3.8")
library(AnnotationHub)
library(GOSemSim)

library(parallel)
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)


#Export packages to each cluster (https://stackoverflow.com/questions/33761123/r-how-can-i-export-methods-provided-by-a-package-to-a-psock-cluster)
clusterEvalQ(cl,c(library(AnnotationHub),library(GOSemSim),library(parallel))) # I tried to pass more than 1 library and it worked


ECK_GO_MF=godata('org.EcK12.eg.db',ont="MF",computeIC = F)
ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)
ECK_GO_CC=godata('org.EcK12.eg.db',ont="CC",computeIC = F)

#Export required variables to each cluster (http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/)
clusterExport(cl=cl, c("ECK_GO_MF","ECK_GO_BP","ECK_GO_CC")) # I tried to pass more than 1 variable and it worked

start.time=Sys.time()

Wang_pairwise_similarity_BP=parLapply(cl=cl,X=pairwise_id_BP,fun=function(annots){
  mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_BP, measure="Wang",combine="BMA")
}
)

Wang_pairwise_similarity_MF=parLapply(cl=cl,X=pairwise_id_MF,fun=function(annots){
  mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_MF, measure="Wang",combine="BMA")
}
) 

Wang_pairwise_similarity_CC=parLapply(cl=cl,X=pairwise_id_CC,fun=function(annots){
  mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_CC, measure="Wang",combine="BMA")
}
) 


save(Wang_pairwise_similarity_BP,Wang_pairwise_similarity_MF,Wang_pairwise_similarity_CC,file="Data/sourced/Wang_pairwise_similarity_separate.RData") #Time difference of 8.254868 mins

end.time=Sys.time()
end.time-start.time 
stopCluster(cl)




strain1_strain2_WangBP_WangMF_WangCC=cbind(strain1_strain2,
                                           BP=unlist(Wang_pairwise_similarity_BP),
                                           MF=unlist(Wang_pairwise_similarity_MF),
                                           CC=unlist(Wang_pairwise_similarity_CC))


##test if they are correct:
##GO1=id_allAttributes[id_allAttributes$ids=="1","GO"] %>% unique
##GO1_BP=GO1[AnnotationDbi::Ontology(GO1)=="BP" & !is.na(AnnotationDbi::Ontology(GO1)=="BP")]

##GO3=id_allAttributes[id_allAttributes$ids=="3","GO"] %>% unique
##GO3_BP=GO3[AnnotationDbi::Ontology(GO3)=="BP" & !is.na(AnnotationDbi::Ontology(GO3)=="BP")]

##mgoSim(GO1_BP,GO3_BP,semData = ECK_GO_BP, measure="Wang",combine="BMA") #0.298


save(strain1_strain2_WangBP_WangMF_WangCC,file="Data/sourced/strain1_strain2_WangBP_WangMF_WangCC.RData")


}








#Distribution of GO pairwise similarity (Wang method - graph based) in Nichols' (Distance = NA is removed)
#Why can a distance be NA? -> a leat one of the strains is NA (it doesn't have GO annotations)
hist(main="BP",
     xlab="Similarity (min=0, max=1)",
     strain1_strain2_WangBP_WangMF_WangCC$BP[
       !is.na(strain1_strain2_WangBP_WangMF_WangCC$BP)
       ])

hist(main="MF",
     xlab="Similarity (min=0, max=1)",
     strain1_strain2_WangBP_WangMF_WangCC$MF[
       !is.na(strain1_strain2_WangBP_WangMF_WangCC$MF)
       ])


hist(main="CC",
     xlab="Similarity (min=0, max=1)",
     strain1_strain2_WangBP_WangMF_WangCC$CC[
       !is.na(strain1_strain2_WangBP_WangMF_WangCC$CC)
       ])







#Correlation VS annotation (here annotations are gene similarity by GO annotations)


#join pcc based distance and sort
strain1strain2_pcc=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","pcc")]
strain1_strain2_WangBP_WangMF_WangCC_pcc=left_join(strain1_strain2_WangBP_WangMF_WangCC,strain1strain2_pcc)

strain1_strain2_WangBP_WangMF_WangCC_pcc=strain1_strain2_WangBP_WangMF_WangCC_pcc[order(
  strain1_strain2_WangBP_WangMF_WangCC_pcc$pcc
),]



## pairs of NA similarity are removed

#strain1_strain2_pcc_WangBP_noNA=strain1_strain2_WangBP_WangMF_WangCC_pcc[,c("strain1","strain2","pcc","BP")]
#strain1_strain2_pcc_WangBP_noNA=strain1_strain2_pcc_WangBP_noNA[!is.na(strain1_strain2_pcc_WangBP_noNA$BP),]
#save(strain1_strain2_pcc_WangBP_noNA,file="Data/sourced/strain1_strain2_pcc_WangBP_noNA.RData")


#strain1_strain2_pcc_WangMF_noNA=strain1_strain2_WangBP_WangMF_WangCC_pcc[,c("strain1","strain2","pcc","MF")]
#strain1_strain2_pcc_WangMF_noNA=strain1_strain2_pcc_WangMF_noNA[!is.na(strain1_strain2_pcc_WangMF_noNA$MF),]
#save(strain1_strain2_pcc_WangMF_noNA,file="Data/sourced/strain1_strain2_pcc_WangMF_noNA.RData")

#strain1_strain2_pcc_WangCC_noNA=strain1_strain2_WangBP_WangMF_WangCC_pcc[,c("strain1","strain2","pcc","CC")]
#strain1_strain2_pcc_WangCC_noNA=strain1_strain2_pcc_WangCC_noNA[!is.na(strain1_strain2_pcc_WangCC_noNA$CC),]
#save(strain1_strain2_pcc_WangCC_noNA,file="Data/sourced/strain1_strain2_pcc_WangCC_noNA.RData")


#quick way to check the result - BP
#=================================================

#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr VS. Annot(GO similarity - BP)',
     xlab="high |PCC| -- rankings of strain pairs  -- low |PCC|",
     ylab="Cumsum of GO similarity (identical=1, complete different=0)",
     (1:dim(strain1_strain2_pcc_WangBP_noNA)[1])[samples],cumsum(strain1_strain2_pcc_WangBP_noNA$BP[samples]),
     type="l")
abline(a=0,b=sum(strain1_strain2_pcc_WangBP_noNA$BP)/dim(strain1_strain2_pcc_WangBP_noNA)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
## => No difference than a randomized order


#=================================================


#quick way to check the result - MF
#=================================================

#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr VS. Annot(GO similarity - MF)',
     xlab="high |PCC| -- rankings of strain pairs  -- low |PCC|",
     ylab="Cumsum of GO similarity (identical=1, complete different=0)",
     (1:dim(strain1_strain2_pcc_WangMF_noNA)[1])[samples],cumsum(strain1_strain2_pcc_WangMF_noNA$MF[samples]),
     type="l")
abline(a=0,b=sum(strain1_strain2_pcc_WangMF_noNA$MF)/dim(strain1_strain2_pcc_WangMF_noNA)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
## => No difference than a randomized order


#=================================================


#quick way to check the result - CC
#=================================================


#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr VS. Annot(GO similarity - CC)',
     xlab="high |PCC| -- rankings of strain pairs  -- low |PCC|",
     ylab="Cumsum of GO similarity (identical=1, complete different=0)",
     (1:dim(strain1_strain2_pcc_WangCC_noNA)[1])[samples],cumsum(strain1_strain2_pcc_WangCC_noNA$CC[samples]),
     type="l")
abline(a=0,b=sum(strain1_strain2_pcc_WangCC_noNA$CC)/dim(strain1_strain2_pcc_WangCC_noNA)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
## => No difference than a randomized order


#=================================================
#Conclusion: More similar profile (by PCC) => higher GO similarity




##Revert the above experiment and ask: Do ranked GO similarity (x-axis) -- cummulative |PCC| (y-axis) + several examples



#GO similarity VS. cumulative |PCC|

#BP
#=================================================
dim(strain1_strain2_pcc_WangBP_noNA)
orderedBy_BP=strain1_strain2_pcc_WangBP_noNA[order(strain1_strain2_pcc_WangBP_noNA$BP,decreasing=T),]

samples=1:4000 #note that this masks samples() function
plot(main='BP',
     xlab="Most similar -- rankings of GO similarity  -- least similar",
     ylab="Cumsum of |PCC|",
     (1:dim(orderedBy_BP)[1])[samples],cumsum(1-orderedBy_BP$pcc[samples]), #in $pcc, it's actually 1-|PCC| => Have to convert it back
     type="l")
abline(a=0,b=sum(1-orderedBy_BP$pcc)/dim(orderedBy_BP)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
#=================================================


#MF
#=================================================
dim(strain1_strain2_pcc_WangMF_noNA)
orderedBy_MF=strain1_strain2_pcc_WangMF_noNA[order(strain1_strain2_pcc_WangMF_noNA$MF,decreasing=T),]

samples=1:4000 #note that this masks samples() function
plot(main='MF',
     xlab="Most similar -- rankings of GO similarity  -- least similar",
     ylab="Cumsum of |PCC|",
     (1:dim(orderedBy_MF)[1])[samples],cumsum(1-orderedBy_MF$pcc[samples]), #in $pcc, it's actually 1-|PCC| => Have to convert it back
     type="l")
abline(a=0,b=sum(1-orderedBy_MF$pcc)/dim(orderedBy_MF)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
#=================================================

#CC
#=================================================
dim(strain1_strain2_pcc_WangCC_noNA)
orderedBy_CC=strain1_strain2_pcc_WangCC_noNA[order(strain1_strain2_pcc_WangCC_noNA$CC,decreasing=T),]

samples=1:4000 #note that this masks samples() function
plot(main='CC',
     xlab="Most similar -- rankings of GO similarity  -- least similar",
     ylab="Cumsum of |PCC|",
     (1:dim(orderedBy_CC)[1])[samples],cumsum(1-orderedBy_CC$pcc[samples]), #in $pcc, it's actually 1-|PCC| => Have to convert it back
     type="l")
abline(a=0,b=sum(1-orderedBy_CC$pcc)/dim(orderedBy_CC)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
#=================================================

#Conclusion: Higher |PCC| -> more similar of 2 sets of GO annotations














#filter GO annotations so only BP terms are left

#retrieve go terms , plot histogram, and maybe wordcloud?



##MF
dim(orderedBy_MF) #3088855*0.05 ~ 154442 (similarity=0.815)

#CC
dim(orderedBy_CC) #2579856*0.05 ~ 128992 (similarity=1)




#If I use 5% and categorize them as co-annotated, others are not, what happens?

#quick way to check the result - BP
#=================================================

strain1_strain2_pcc_binaryWangBP_noNA=strain1_strain2_pcc_WangBP_noNA
strain1_strain2_pcc_binaryWangBP_noNA$BP=ifelse(strain1_strain2_pcc_binaryWangBP_noNA$BP>=orderedBy_BP$BP[173383],1,0)


#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr VS. Annot(BP)',
     xlab="high |PCC| -- rankings of strain pairs  -- low |PCC|",
     ylab="Cumsum of co-annotation",
     (1:dim(strain1_strain2_pcc_binaryWangBP_noNA)[1])[samples],cumsum(strain1_strain2_pcc_binaryWangBP_noNA$BP[samples]),
     type="l")
abline(a=0,b=sum(strain1_strain2_pcc_binaryWangBP_noNA$BP)/dim(strain1_strain2_pcc_binaryWangBP_noNA)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
## => No difference than a randomized order


#=================================================


#quick way to check the result - MF
#=================================================
strain1_strain2_pcc_binaryWangMF_noNA=strain1_strain2_pcc_WangMF_noNA
strain1_strain2_pcc_binaryWangMF_noNA$MF=ifelse(strain1_strain2_pcc_binaryWangMF_noNA$MF>=orderedBy_MF$MF[154442],1,0)


#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr VS. Annot(MF)',
     xlab="high |PCC| -- rankings of strain pairs  -- low |PCC|",
     ylab="Cumsum of co-annotation",
     (1:dim(strain1_strain2_pcc_binaryWangMF_noNA)[1])[samples],cumsum(strain1_strain2_pcc_binaryWangMF_noNA$MF[samples]),
     type="l")
abline(a=0,b=sum(strain1_strain2_pcc_binaryWangMF_noNA$MF)/dim(strain1_strain2_pcc_binaryWangMF_noNA)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
## => No difference than a randomized order


#=================================================


#quick way to check the result - CC
#=================================================
## pairs of NA similarity are removed
strain1_strain2_pcc_binaryWangCC_noNA=strain1_strain2_pcc_WangCC_noNA
strain1_strain2_pcc_binaryWangCC_noNA$CC=ifelse(strain1_strain2_pcc_binaryWangCC_noNA$CC>=orderedBy_CC$CC[128992],1,0)


#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr VS. Annot(CC)',
     xlab="high |PCC| -- rankings of strain pairs  -- low |PCC|",
     ylab="Cumsum of co-annotation",
     (1:dim(strain1_strain2_pcc_binaryWangCC_noNA)[1])[samples],cumsum(strain1_strain2_pcc_binaryWangCC_noNA$CC[samples]),
     type="l")
abline(a=0,b=sum(strain1_strain2_pcc_binaryWangCC_noNA$CC)/dim(strain1_strain2_pcc_binaryWangCC_noNA)[1],col="red",lty="dotted") #random expectation (avg |PCC|)
## => No difference than a randomized order

#=================================================
#Conclusion: using 5% as the cutoff for co-annotation have similar result compared to using cumulative GO similarity







