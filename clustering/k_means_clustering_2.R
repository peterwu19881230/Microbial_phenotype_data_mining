#Ref: https://www.datacamp.com/community/tutorials/k-means-clustering-r

#Perform k-means. If the current seed has a convergence issue => seed=seed+1 => redo until there's no convergence issue
safe_k_means=function(df,k,seed){
  set.seed(seed)
  tryCatch(kmeans(df,k)$cluster, warning=function(war){ ##Ref: http://mazamascience.com/WorkingWithData/?p=912
    seed=seed+1 #this looks for the seed argument in the above level
    safe_k_means(df,k,seed)
  })
}

##code to test the function safe_k_means:
##set.seed(102); kmeans(All_Data_NAimputed,17) #Warning message: did not converge in 10 iterations 
##temp=safe_k_means(All_Data_NAimputed,17,102)
##set.seed(103); temp2=kmeans(All_Data_NAimputed,17)$cluster #Warning message: did not converge in 10 iterations 
##set.seed(104); temp2=kmeans(All_Data_NAimputed,17)$cluster 
##identical(temp,temp2)


cluster_1to50_df=data.frame(id=1:3979,k_1=rep(1,3979)) #add a column for k=1
for(i in 2:50){
  cluster_1to50_df=cbind(cluster_1to50_df,safe_k_means(df=All_Data_NAimputed,k=i,seed=102))
}


names(cluster_1to50_df)=c("id",paste("k_",rep(1:50),sep=""))





write.table(cluster_1to50_df,file="Data/k_means_Nichols_1to50.txt",quote=F,sep="\t",row.names = F)


