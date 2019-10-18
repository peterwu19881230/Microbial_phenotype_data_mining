#Randomize the strain IDs + Use all annotated terms to genes to validate k-means clustering

#Calculate pairwise semantic similarity of GO BP terms in k means
#1. all pairwise similarity for all unique terms for each cluster in each k
#2. all pairwise similarity for all unique terms for all genes
#3. analysis




library(foreach)
library(doParallel)
library(parallel)
library(doSNOW)
library(pbapply)

#goea=read.table("Data/k_means_Nichols_1to50.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
goea=read.table("Data/k_means_Nichols.txt",header=T,sep="\t",quote="",stringsAsFactors=F)






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


result_list=lapply(6:12,function(Col){#iterate through all k value,starting from the 2nd column (for this script I do 6~12) 
  
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
    
    if(length(strainID)>=2){
      strainComb=t(combn(strainID,2))
      
      sim=c()
      for(i in 1:dim(strainComb)[1]){ #Calculate pairwise similarity for all the terms
        
        go1=id_allGO$GO[id_allGO$ids %in% strainComb[i,1]] #all terms for gene 1
        go2=id_allGO$GO[id_allGO$ids %in% strainComb[i,2]] #all terms for gene 2
        
        
        sim[i]=GOSemSim::mgoSim(go1, go2, semData=ECK_GO_BP,measure="Wang",combine="BMA")
        
        
        
        simDF=rbind(simDF,data.frame(strain1=strainComb[i,1],
                                     strain2=strainComb[i,2],
                                     cluster_id=cluster,
                                     k=k_val[Col-1],
                                     sim=sim[i],
                                     avgSim=NA, #This is the average similarity within that k (value will be appended later)
                                     medianSim=NA #This is the median similarity within that k (value will be appended later)
        ))
        
        
      }
      
      simDF$avgSim=mean(sim,na.rm=T)
      simDF$medianSim=median(sim,na.rm=T)
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



GOEA_semantic_DF_geneSim=resultDF
#save(GOEA_semantic_DF_geneSim,file="Data/GOEA_semantic_DF_geneSim.RData")


#2. all pairwise similarity for all genes
allBPMean=mean(strain1strain2_allAnnotations_allDistances$Wang_BP,na.rm=T)




#3. analysis

##average of similarity of clusters for each k 

##this calculates means based on equal weight of each cluster
'
Table=GOEA_semantic_DF_geneSim[,c("cluster_id","k","avgSim")] %>% unique

cluster_mean_forEachK=c()
for(i in 1:length(unique(Table$k))){

k=unique(Table$k)[i]
cluster_mean_forEachK[i]=mean(Table$avgSim[Table$k==k],na.rm=T)
}

names(cluster_mean_forEachK)=unique(Table$k)
'


##this calculates means based on sizes of clusters
Table=GOEA_semantic_DF_geneSim[,c("k","sim")] 

cluster_mean_forEachK=c()
for(i in 1:length(unique(Table$k))){
  k=Table$k[i]
  cluster_mean_forEachK[i]=mean(Table$sim[Table$k==k],na.rm=T)
}

names(cluster_mean_forEachK)=unique(Table$k)



#plot
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
pdf(file=paste(currentDir,"kmeans_GOBPsemantic_geneSim.pdf",sep="/"),width=8, height=5)

barplot(cluster_mean_forEachK,ylim=c(0,1),xlab = "k in k-means",ylab="GO BP similarity") #average of average similarity within each cluster for each k
abline(a=allBPMean,b=0,lty="dotted") #average of all pairwise similarity of BP terms


dev.off()


