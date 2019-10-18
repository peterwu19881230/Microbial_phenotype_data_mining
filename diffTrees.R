#Goal: Determine whether the all.equal(,use.edge.length = FALSE) can tell the difference when the topology of the trees are different


#Clustering of original data
mat.fig.4A.strain<-All_Data[which(ECK_1st_table$sorted_ECK_missing_gene_names_added %in% c("ECK0457-ACRA","ECK0456-ACRB","ECK1524-MARA","ECK1525-MARB","ECK1523-MARR")),]
dist<-dist(mat.fig.4A.strain,method="euclidean")
hclust_original<-hclust(dist,method="complete")
dendro5genes_quantitative<-as.dendrogram(hclust_original)

#Clustering of the binary data: 5 genotypes from fig.4A using all 324 conditions 
mat.binary.fig.4A.strain<-Binary_Data[which(ECK_1st_table$sorted_ECK_missing_gene_names_added %in% c("ECK0457-ACRA","ECK0456-ACRB","ECK1524-MARA","ECK1525-MARB","ECK1523-MARR")),]
dist<-dist(mat.binary.fig.4A.strain,method="euclidean")
hclust_5genes<-hclust(dist,method="complete")
dendro5genes<-as.dendrogram(hclust_5genes)

#dput(dendro5genes,"dendro5genes.txt")   => See the details inside the dendro in order to change the order of it

#Different dendros: 



#plot them
par(mfrow = c(1,2))
dendro5genes %>% plot
dendro5genes_quantitative %>% plot








