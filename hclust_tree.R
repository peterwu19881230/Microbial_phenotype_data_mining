#Generate hierarchical clustering follwing: http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm  by using Nichols' data

#Import the data from .csv and tell the read.csv() that 1st column, 1st row are there for the names of the data
Data<-read.csv("allData_rows_fixed.csv",row.names = 1,header = TRUE)


#Transpose for cor() to work
t_Data<-t(Data)


#Do clustering
withoutNAs_cor <- cor(t_Data,use="pairwise.complete.obs")
dissimilarity <- 1 - withoutNAs_cor
distance <- as.dist(dissimilarity)

pdf('hclustering_genotype')
hclust(distance)
