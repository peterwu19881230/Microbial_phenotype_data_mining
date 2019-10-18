#Calculate pairwise semantic similarity of GO BP terms in k means
#1. all pairwise similarity for all enriched terms for each cluster
#2. pairwise similarity for enriched terms within clusters for each k
#3. analysis


library(GOSemSim)
library(org.EcK12.eg.db)


#goea=read.csv("Data/goea_1to50.csv")
goea=read.table("Data/goea_k5to3500.txt",header=T,sep="\t",quote="",stringsAsFactors=F)



##1. and 2. done together
k_val=goea$k_val %>% unique
ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)

start.time=Sys.time()

simDF=data.frame()
sim_in_clusterDF=data.frame()
for(k in k_val){ #iterate through all k value  
  
  print(paste("Now doing k= ",k,sep=""))
  
  GOterms=goea$GO[goea$NS=="BP" & goea$NS!="N/A" #this filters terms that are NA
                  & goea$k_val==k]
  
  if(length(GOterms)>=2){ #There has to be at least 2 GO terms to do pairwise similarity
    
    comb=t(combn(length(GOterms),2))
    sim=c()
    for(i in 1:dim(comb)[1]){ #Calculate pairwise similarity for all the terms
      
      
      ##1. pairwise similarity for all combination of terms
      go1=GOterms[comb[i,1]]
      go2=GOterms[comb[i,2]]
      
      sim[i]=goSim(go1, go2, semData=ECK_GO_BP,measure="Wang")
      
      simDF=rbind(simDF,data.frame(go1=go1,go2=go2,
                                   k=k,
                                   sim=sim[i],
                                   avgSim=NA, #This is the average similarity within that k (value will be appended later)
                                   medianSim=NA #This is the median similarity within that k (value will be appended later)
      ))
    }
    simDF$avgSim[simDF$k==k]=mean(sim)
    simDF$medianSim[simDF$k==k]=median(sim)
    
  } 
  
  ##2. pairwise similarity for terms within each cluster
  cluster_ids=goea$cluster_id[goea$k_val==k] %>% unique
  
  
  
  for(cluster_id in cluster_ids){
    GOterms=goea$GO[goea$NS=="BP" & goea$NS!="N/A" #this filters terms that are NA
                    & goea$k_val==k 
                    & goea$cluster_id==cluster_id]
    
    if(length(GOterms)>=2){ #There has to be at least 2 GO terms to do pairwise similarity
      
      comb=t(combn(length(GOterms),2))
      sim=c()
      
      for(i in 1:dim(comb)[1]){ #Calculate pairwise similarity for all the terms
        
        
        ##1. pairwise similarity for all combination of terms
        go1=GOterms[comb[i,1]]
        go2=GOterms[comb[i,2]]
        
        sim[i]=goSim(go1, go2, semData=ECK_GO_BP)
        
        sim_in_clusterDF=rbind(sim_in_clusterDF,data.frame(go1=go1,go2=go2,
                                                           k=k,
                                                           cluster_id=cluster_id,
                                                           sim=sim[i],
                                                           avgSim=NA, #This is the average similarity within that k (value will be appended later)
                                                           medianSim=NA #This is the median similarity within that k (value will be appended later)
        ))
      }
      sim_in_clusterDF$avgSim[sim_in_clusterDF$k==k & sim_in_clusterDF$cluster_id==cluster_id]=mean(sim)
      sim_in_clusterDF$medianSim[sim_in_clusterDF$k==k & sim_in_clusterDF$cluster_id==cluster_id]=median(sim)
      
      
    }
    
  }
  
}

end.time=Sys.time()
end.time-start.time 


#average similarity for each k
mean_for_each_k=c()
i=1
for(k in unique(simDF$k)){
  mean_for_each_k[i]=simDF$avgSim[simDF$k==k] %>% unique 
  i=i+1
}


#average of average similarity within each cluster for each k
cluster_mean_forEachK=c()
for(k in unique(sim_in_clusterDF$k)){
  means_for_k=c()
  for(cluster_id in unique(sim_in_clusterDF$cluster_id[sim_in_clusterDF$k==k])){
    means_for_k=c(means_for_k,
                  sim_in_clusterDF$avgSim[sim_in_clusterDF$k==k & sim_in_clusterDF$cluster_id==cluster_id] %>% unique )
    
    
    
  }
  cluster_mean_forEachK=c(cluster_mean_forEachK,mean(means_for_k))
}





#3. analysis
#stacked barplot to show the difference 

blue = rgb(0, 0, 1, alpha=0.2)
red = rgb(1, 0, 0, alpha=0.2)

barplot(cluster_mean_forEachK,ylim=c(0,1),col=red,xlab = "k in k-means",ylab="GO BP similarity") #average of average similarity within each cluster for each k
barplot(mean_for_each_k,ylim=c(0,1),col=blue,add=T) #average similarity for each k


#This script takes more than 1 day on my PC and still hasn't finished. Need optimization. 
#However, I feel reluctantly to do so because the result doesn't look like a fig we can put on a paper 
#(From the partial result (K<=45) I got after terminating the process, average of average similarity within each cluster for each k
#is always greater than average similarity for each k




