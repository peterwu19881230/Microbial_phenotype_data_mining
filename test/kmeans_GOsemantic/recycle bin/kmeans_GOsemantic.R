#Calculate pairwise semantic similarity of GO BP terms in k means
#1. all pairwise similarity for all enriched terms for each cluster
#2. pairwise similarity for enriched terms within clusters for each k
#3. analysis


library(GOSemSim)
library(org.EcK12.eg.db)

library(foreach)
library(doParallel)
library(parallel)

#goea=read.csv("Data/goea_1to50.csv")
goea=read.table("Data/goea_k5to3500.txt",header=T,sep="\t",quote="",stringsAsFactors=F)



##1. and 2. done together
k_val=goea$k_val %>% unique
ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)

start.time=Sys.time()

numCores=detectCores()-1
cl=makeCluster(numCores)


result_list=foreach(k = k_val) %dopar% { #iterate through all k value  
  
  simDF=data.frame()
  
  print(paste("Now doing k= ",k,sep=""))
  
  GOterms=goea$GO[goea$NS=="BP" & goea$NS!="N/A" #this filters terms that are NA
                  & goea$k_val==k]
  
  cluster_ids=goea$cluster_id[goea$NS=="BP" & goea$NS!="N/A" #this filters terms that are NA
                              & goea$k_val==k]
  
  
  if(length(GOterms)>=2){ #There has to be at least 2 GO terms to do pairwise similarity
    
    comb=t(combn(length(GOterms),2))
    sim=c()
    for(i in 1:dim(comb)[1]){ #Calculate pairwise similarity for all the terms
      
      
      
      go1=GOterms[comb[i,1]]; cluster_id1=cluster_ids[comb[i,1]]
      go2=GOterms[comb[i,2]]; cluster_id2=cluster_ids[comb[i,2]]
      
      
      sim[i]=goSim(go1, go2, semData=ECK_GO_BP,measure="Wang")
      
      simDF=rbind(simDF,data.frame(go1=go1,cluster_id1=cluster_id1,
                                   go2=go2,cluster_id2=cluster_id2,
                                   k=k,
                                   sim=sim[i],
                                   avgSim=NA, #This is the average similarity within that k (value will be appended later)
                                   medianSim=NA #This is the median similarity within that k (value will be appended later)
      ))
    }
    simDF$avgSim[simDF$k==k]=mean(sim)
    simDF$medianSim[simDF$k==k]=median(sim)
    
  } 
  
  
  simDF
}

end.time=Sys.time()
end.time-start.time 


##If the result is empty for that k, remove it
index_for_NULL=c()
for(i in 1:length(result_list)){
  if(identical(result_list[[i]],data.frame())){ #(if the result is empty for that k)
    index_for_NULL=c(index_for_NULL,i)
  }
}


new_result_list=result_list[-index_for_NULL]


##Bind all dataframe results into one
resultDF=data.frame()
for(i in 1:length(new_result_list)){
  resultDF=rbind(resultDF,new_result_list[[i]])
}


GOEA_semantic_DF=resultDF

#save(GOEA_semantic_DF,file="Data/GOEA_semantic_DF")
load("Data/GOEA_semantic_DF")


#average similarity for each k
mean_for_each_k=c()
i=1
for(k in unique(GOEA_semantic_DF$k)){
  mean_for_each_k[i]=GOEA_semantic_DF$avgSim[GOEA_semantic_DF$k==k] %>% unique 
  i=i+1
}

names(mean_for_each_k)=unique(GOEA_semantic_DF$k) #give each mean the k where it was originally from


#average of average similarity within each cluster for each k
cluster_mean_forEachK=c()
for(k in unique(GOEA_semantic_DF$k)){
  
  TFvec=ifelse(GOEA_semantic_DF$k==k & GOEA_semantic_DF$cluster_id1==GOEA_semantic_DF$cluster_id2,T,F)
  
  inClusterDF=GOEA_semantic_DF[TFvec,]
  
  means_for_k=c()
  for(cluster_id in unique(inClusterDF$cluster_id1)){ #Using cluster_id1 is enough. Don't have to use cluster_id2
    means_for_k=c(means_for_k,
                  inClusterDF[inClusterDF$cluster_id1==cluster_id,]$sim %>% mean)
  }
  
  cluster_mean_forEachK=c(cluster_mean_forEachK,mean(means_for_k))
}


names(cluster_mean_forEachK)=unique(GOEA_semantic_DF$k) #give each mean the k where it was originally from




#3. analysis
#stacked barplot to show the difference 

currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
pdf(file=paste(currentDir,"kmeans_GOBPsemantic.pdf",sep="/"),width=8, height=5)

blue = rgb(0, 0, 1, alpha=0.2)
red = rgb(1, 0, 0, alpha=0.2)

barplot(cluster_mean_forEachK,ylim=c(0,1),col=red,xlab = "k in k-means",ylab="GO BP similarity") #average of average similarity within each cluster for each k
barplot(mean_for_each_k,ylim=c(0,1),col=blue,add=T) #average similarity for each k

dev.off()


