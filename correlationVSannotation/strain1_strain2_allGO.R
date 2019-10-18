#Pairwise semantic similarity of Nichols strains (by sets of GO annotations (BP, MF, CC all included) )


#How many ids in Nichols got GO annotations?
names(id_allAttributes)
id_GO=id_allAttributes[,c("ids","GO")] %>% unique
unique(id_GO$ids[!is.na(id_GO$GO)]) %>% length #3306



#BiocManager::install("AnnotationHub", version = "3.8")
library(AnnotationHub)
library(GOSemSim)

#1st and 2nd column for strain pairs
##strain1_strain2=as.data.frame(t(combn(1:3979,2)))
##names(strain1_strain2)=c("strain1","strain2")
##strain1_strain2$strain1=as.character(strain1_strain2$strain1)
##strain1_strain2$strain2=as.character(strain1_strain2$strain2)
##save(strain1_strain2,file="Data/sourced/strain1_strain2.RData")



if(!file.exists("Data/sourced/GOannotation.RData")){
  strain1strain2_transposed=as.data.frame((combn(1:3979,2))) 
  
  
  start.time=Sys.time()
  
  library(parallel)
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  
  clusterExport(cl=cl,"id_GO") # I tried to pass more than 1 variable and it worked
  
  GOannotation=parLapply(cl=cl,strain1strain2_transposed,fun=function(strain1strain2){ 
    ##Although input is not a list, it will be coerced by as.list() according to ?lapply
    
    strain1=as.character(strain1strain2[1])
    strain2=as.character(strain1strain2[2])
    
    annotation1=id_GO$GO[id_GO$ids==strain1]
    annotation2=id_GO$GO[id_GO$ids==strain2]
    
    return(list(annotation1,annotation2))
  })
  
  end.time=Sys.time()
  end.time-start.time #Time difference of 11.35332 mins on my PC if done by 7 cpus
  
  save(GOannotation,file="Data/sourced/GOannotation.RData")
  stopCluster(cl)
}


str(GOannotation)




#Get combined (BP, MF, CC altogether) distance 
if(!file.exists("Data/Wang_pairwise_similarity.RData")){
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
  
  Wang_pairwise_similarity=parLapply(cl=cl,X=GOannotation,fun=function(annots){
    MF=mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_MF, measure="Wang",combine="BMA")
    BP=mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_BP, measure="Wang",combine="BMA")
    CC=mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_CC, measure="Wang",combine="BMA")
    return(list(MF=MF,BP=BP,CC=CC)) 
    
  }
  ) 
  
  save(Wang_pairwise_similarity,file="Data/sourced/Wang_pairwise_similarity.RData")
  
  end.time=Sys.time()
  end.time-start.time #Time difference of 43.23664 mins from my PC if done by 7 cpus
  stopCluster(cl)
  
  #==non-parallelized code example==
  #ECK_GO_MF=godata('org.EcK12.eg.db',ont="MF",computeIC = F)
  #ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)
  #ECK_GO_CC=godata('org.EcK12.eg.db',ont="CC",computeIC = F)
  #Wang_distance=sapply(X=GOannotation[1:50],FUN=function(annots){
  #  MF=mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_MF, measure="Wang",combine="BMA")
  #  BP=mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_BP, measure="Wang",combine="BMA")
  #  CC=mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_CC, measure="Wang",combine="BMA")
  #  return(list(MF=MF,BP=BP,CC=CC)) 
  #}
  #) 
  
  #Probelms:
  #MF,BP,CC seem to be giving all the same results.  
  #For the GO annotations for Nichols' strains, they are mixture of all MF, BP ,CC 
  #=> Does it make sense that I use all of the annotations to get combinatorial results for each strain pairs?
}



##Did MF, BP, CC all give me the same results?=> Yes
'
start.time=Sys.time()

simplified_Wang_pairwise_distance_similarity=sapply(Wang_pairwise_similarity,FUN = function(similarity){
unique(unlist(similarity)) %>% length   
})

end.time=Sys.time()
end.time-start.time #Time difference of 6.949467 mins

sum(simplified_Wang_pairwise_similarity_length==1) #7914231 => all of the gene pairs have same distance regardless of using MF, BP, CC semantic data
'




## Simplify Wang distance since MF, BP, CC all give me the same results => Giving 1 value is enough
start.time=Sys.time()

simplified_Wang_pairwise_similarity=sapply(Wang_pairwise_similarity,FUN = function(similarity){
  unique(unlist(similarity))  
})

end.time=Sys.time()
end.time-start.time #1.455671 mins




strain1_strain2_WangSimilarity=cbind(strain1_strain2,simplified_Wang_pairwise_similarity)
head(strain1_strain2_WangSimilarity)

summary(strain1_strain2_WangSimilarity$simplified_Wang_pairwise_similarity) #Why are there so many NA values (2451066)?

##No. of similarity values caused by no annotation in at least 1 strain: 1 NA (no GO annotation to that strain) + 2 NAs
no_oof_NA=3979-3306
(3979-no_oof_NA)*no_oof_NA+choose(no_oof_NA,2) #2451066 (I have double checked that this is the correct combination by using a smaller vector)
##The number (2451066) matches => All NAs from Wang distance(similarity) resulted from at least 1 strain having no GO annotations

hist(strain1_strain2_WangSimilarity$simplified_Wang_pairwise_similarity[
  !is.na(strain1_strain2_WangSimilarity$simplified_Wang_pairwise_similarity)
  ])

#Correlation VS annotation (here annotations are gene similarity by GO annotations)

## Remove NA from strain1_strain2_WangSimilarity$simplified_Wang_pairwise_similarity
## pairs of NA similarity are removed
strain1_strain2_WangSimilarity_noNA=strain1_strain2_WangSimilarity[!is.na(strain1_strain2_WangSimilarity$simplified_Wang_pairwise_similarity),]

names(strain1strain2_allAnnotations_allDistances)
strain1strain2_pcc=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","pcc")]


#join pcc based distance and sort
strain1_strain2_WangSimilarity_noNA_pcc=left_join(strain1_strain2_WangSimilarity_noNA,strain1strain2_pcc)
strain1_strain2_WangSimilarity_noNA_pcc=strain1_strain2_WangSimilarity_noNA_pcc[order(
  strain1_strain2_WangSimilarity_noNA_pcc$pcc
),]


#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr VS. Annot(GO similarity - MF, BP, CC combined)',
     xlab="high |PCC| -- rankings of strain pairs  -- low |PCC|",
     ylab="Cumsum of GO similarity (identical=1, complete different=0)",
     (1:dim(strain1_strain2_WangSimilarity_noNA_pcc)[1])[samples],cumsum(strain1_strain2_WangSimilarity_noNA_pcc$simplified_Wang_pairwise_similarity[samples]),
     type="l")
## => A little difference than a randomized order. Not sure how to calculate the significance








