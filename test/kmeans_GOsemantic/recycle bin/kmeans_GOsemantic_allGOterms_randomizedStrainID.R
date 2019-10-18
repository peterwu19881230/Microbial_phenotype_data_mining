#Randomize the strain IDs + Use all annotated terms to genes to validate k-means clustering

#Calculate pairwise semantic similarity of GO BP terms in k means
#1. all pairwise similarity for all unique terms for each cluster in each k
#2. all pairwise similarity for all unique terms for all genes
#3. analysis




library(foreach)
library(doParallel)
library(parallel)
library(doSNOW)


#goea=read.table("Data/k_means_Nichols_1to50.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
goea=read.table("Data/k_means_Nichols.txt",header=T,sep="\t",quote="",stringsAsFactors=F)


#scramble the ids and clusters of strains
set.seed(102)
for(i in 1:dim(goea)[2]){
  goea[,i]=sample(goea[,i])
}



id_allGO=id_allAttributes[,c("ids","GO")] %>% unique
id_allGO=id_allGO[!is.na(id_allGO$GO),]
id_allGO=id_allGO[AnnotationDbi::Ontology(id_allGO$GO)=="BP" & !is.na(AnnotationDbi::Ontology(id_allGO$GO)=="BP"),]



#1. all pairwise similarity for all unique terms for each cluster in each k

library(GOSemSim)
library(org.EcK12.eg.db)
k_val=c(5,10,20,50,100,200,500,1000,2000,3000,3500) 
## for k=5,10,20,50 (especially k=5) it takes really long time (probably more than 1 day).
## I am not going to run it because from the result of the rest ks I already know the trend

ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)


start.time=Sys.time()


result_list=lapply(11:12,function(Col){#iterate through all k value,starting from the 2nd column (for this script I do 6~12) 
  
  numCores=detectCores()-1
  cl=makeCluster(numCores)
  registerDoSNOW(cl) #if I don't do this on my PC (windows), things won't be executed in parallel
  
  clusterExport(cl,c("goea","id_allGO","k_val","ECK_GO_BP")) 
  ##I shouldn't export "Col" (because it's already available within the function) and I can't (because clusterExport looks variables in global environment, not within the function)
  
  
  print(paste("Now doing k= ",k_val[Col-1],sep=""))
  
  clusters=1:k_val[Col-1]
  
  simDF_list=pblapply(clusters,function(cluster){
    
    simDF=data.frame()
    
    print(paste("Now doing cluster= ",cluster,sep=""))
    
    #strain IDs in that cluster
    strainID=goea$id[goea[[Col]]==cluster]
    
    #All unique GO terms within that cluster
    uniqueGOterms=unique(id_allGO$GO[id_allGO$ids %in% strainID])
    
    if(length(uniqueGOterms)>=2){ #There has to be at least 2 GO terms to do pairwise similarity
      
      comb=t(combn(length(uniqueGOterms),2))
      sim=c()
      for(i in 1:dim(comb)[1]){ #Calculate pairwise similarity for all the terms
        
        go1=uniqueGOterms[comb[i,1]]
        go2=uniqueGOterms[comb[i,2]]
        
        sim[i]=GOSemSim::goSim(go1, go2, semData=ECK_GO_BP,measure="Wang")
        
        
        
        simDF=rbind(simDF,data.frame(go1=go1,
                                     go2=go2,
                                     cluster_id=cluster,
                                     k=k_val[Col-1],
                                     sim=sim[i],
                                     avgSim=NA, #This is the average similarity within that k (value will be appended later)
                                     medianSim=NA #This is the median similarity within that k (value will be appended later)
        ))
        
        
      }
      simDF$avgSim=mean(sim)
      simDF$medianSim=median(sim)
      
    }
    
    return(simDF)
    
  },cl=cl)
  
  
  stopCluster(cl) #Why is this not stopping the cpus properly?
  
  simDF=data.frame()
  for(i in 1:length(simDF_list)){
    simDF=rbind(simDF,simDF_list[[i]])
  }
  
  
  
  
  
  return(simDF)
  
}) 






end.time=Sys.time()
end.time-start.time #Time difference of 19.78207 mins


#Bind all dataframe results into one
resultDF=data.frame()
for(i in 1:length(result_list)){
  resultDF=rbind(resultDF,result_list[[i]])  
}



GOEA_semantic_DF_allTerms=resultDF
#save(GOEA_semantic_DF_allTerms,file="Data/GOEA_semantic_DF_allTerms.RData")


#2. all pairwise similarity for all unique terms for all genes
load("Data/GOBP1_GOBP2_similarity.RData")
allBPMean=mean(GOBP1_GOBP2_similarity$similarity)




#3. analysis

##average of similarity of clusters for each k 

##this calculates means based on equal weight of each cluster
'
Table=GOEA_semantic_DF_allTerms[,c("cluster_id","k","avgSim")] %>% unique

cluster_mean_forEachK=c()
for(i in 1:length(unique(Table$k))){
  
  k=unique(Table$k)[i]
  cluster_mean_forEachK[i]=mean(Table$avgSim[Table$k==k])
}

names(cluster_mean_forEachK)=unique(Table$k)
'


##this calculates means based on sizes of clusters
Table=GOEA_semantic_DF_allTerms[,c("k","sim")] 

cluster_mean_forEachK=c()
for(i in 1:length(unique(Table$k))){
  k=Table$k[i]
  cluster_mean_forEachK[i]=mean(Table$sim[Table$k==k])
}

names(cluster_mean_forEachK)=unique(Table$k)



#plot
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
pdf(file=paste(currentDir,"kmeans_GOBPsemantic_all.pdf",sep="/"),width=8, height=5)

barplot(cluster_mean_forEachK,ylim=c(0,1),xlab = "k in k-means",ylab="GO BP similarity") #average of average similarity within each cluster for each k
abline(a=allBPMean,b=0,lty="dotted") #average of all pairwise similarity of BP terms


dev.off()


