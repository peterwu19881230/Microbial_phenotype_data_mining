#Run k=50 multiple times to see how stable k-means is 


#Ref: https://www.datacamp.com/community/tutorials/k-means-clustering-r

start.time=Sys.time()

result=data.frame(1:3979)
seed=1
while(dim(result)[2]!=101){ #Set no. of iternations to 100 (set 101 because the 1st column is the strain name)
  
  set.seed(seed)
  clusters_50=tryCatch({kmeans(All_Data_NAimputed,200)},warning = function(w) { print("warning") })
  
  if(!identical(clusters_50,"warning")){
    result=cbind(result,clusters_50$cluster)
  }
  
  seed=seed+1
  print(seed)
}


colnames(result)=c("id",paste('50 cluster',1:100,sep="_"))

end.time=Sys.time()
end.time-start.time #Time difference of 10.34227 mins



write.table(result,file="Data/k_means_Nichols_50For100Times.txt",quote=F,sep="\t",row.names = F)








