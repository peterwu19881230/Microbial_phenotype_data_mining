#Ref: https://www.datacamp.com/community/tutorials/k-means-clustering-r
set.seed(102)
clusters_5=kmeans(All_Data_NAimputed,5)
set.seed(102)
clusters_10=kmeans(All_Data_NAimputed,10)
set.seed(102)
clusters_20=kmeans(All_Data_NAimputed,20)
set.seed(102)
clusters_50=kmeans(All_Data_NAimputed,50)
set.seed(102)
clusters_100=kmeans(All_Data_NAimputed,100)
set.seed(103)
clusters_200=kmeans(All_Data_NAimputed,200) #showed "did not converge in 10 iterations" when seed=102 so I used a different seed 
set.seed(103)
clusters_500=kmeans(All_Data_NAimputed,500) #showed "did not converge in 10 iterations" when seed=102 so I used a different seed 
set.seed(102)
clusters_1000=kmeans(All_Data_NAimputed,1000)
set.seed(102)
clusters_2000=kmeans(All_Data_NAimputed,2000)
set.seed(102)
clusters_3000=kmeans(All_Data_NAimputed,3000)
set.seed(102)
clusters_3500=kmeans(All_Data_NAimputed,3500)



head(clusters_5$cluster)
length(clusters_5$cluster)


result=rbind(1:3979,
             clusters_5$cluster,
             clusters_10$cluster,
             clusters_20$cluster,
             clusters_50$cluster,
             clusters_100$cluster,
             clusters_200$cluster,
             clusters_500$cluster,
             clusters_1000$cluster,
             clusters_2000$cluster,
             clusters_3000$cluster,
             clusters_3500$cluster
             ) %>% t

colnames(result)=c("id","5 cluster","10 cluster","20 cluster",
                "50 cluster","100 cluster","200 cluster","500 cluster",
                "1000 cluster","2000 cluster","3000 cluster","3500 cluster")
head(result)


write.table(result,file="Data/k_means_Nichols.txt",quote=F,sep="\t",row.names = F)


#Distribution of the no. of genes in each cluster
result=as.data.frame(result)

table(result$'5 cluster')
table(result$'10 cluster')
table(result$'20 cluster')
table(result$'50 cluster')
table(result$'100 cluster')
table(result$'200 cluster')
table(result$'500 cluster')
table(result$'1000 cluster')
table(result$'2000 cluster')
table(result$'3000 cluster')
table(result$'3500 cluster')



