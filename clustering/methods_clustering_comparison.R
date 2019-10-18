#Different methods for comparing clusterings:

library("dendextend")
#FM-index

#The following takes too much time on my mac
#Using euclidean+pcc for both types of data
result<-c()
for(i in 2:3978){ #This number range is all the tree numbers allowed in FM-index when there are 3979 samples
result<-c(result,FM_index(cutree(hclust_euclidean_complete,k=i),cutree(binary_hclust_euclidean_complete,k=i)))
}
FM_index_original_binary_euclidean_complete<-result
rm(result)
#When k=3: FM=index= 0.9624
#When k=10: FM=index= 0.9684
#When k=100: FM=index= 0.685
#When k=1000: FM=index= 0.0861
#When k=3978: FM=index= 0

#Using pcc->dist+complete for original and euclidean+complete for binary
FM_index(cutree(hclust_pcc_complete,k=3),cutree(binary_hclust_euclidean_complete,k=3))
FM_index(cutree(hclust_pcc_complete,k=10),cutree(binary_hclust_euclidean_complete,k=10))
FM_index(cutree(hclust_pcc_complete,k=100),cutree(binary_hclust_euclidean_complete,k=100))
FM_index(cutree(hclust_pcc_complete,k=1000),cutree(binary_hclust_euclidean_complete,k=1000))
FM_index(cutree(hclust_pcc_complete,k=3978),cutree(binary_hclust_euclidean_complete,k=3978))

#When k=3: FM=index= 0.5842812
#When k=10: FM=index= 0.308538
#When k=100: FM=index= 0.1078998
#When k=1000: FM=index= 0.0194877
#When k=3978: FM=index= 0


#Conclusion, if using both euclidean+complete, they look similar only if No. of clusterings are small. 
#For: pcc->dist+complete for original and euclidean+complete for binary => They don't look that similar even if the No. of clustering is small