#Goal: See if euclidean distance can help sort the relavence of ECKs

#Sort the output of a dist obj: https://stackoverflow.com/questions/31591546/sorting-the-output-of-dist
m1=as.matrix(binary_dist_euclidean)
m1[upper.tri(m1)]=0 
library(reshape2)
m2=subset(melt(m1), value!=0)
m2=m2[,c(2,1,3)] ##reorder the columns
names(m2)=c("strain1","strain2","Euclidean Distance") ##name the columns
m2=m2[order(m2[,3]),] ##reorder by Euclidean Distance
sort.binary_dist_euclidean=m2; rm(m2)


#example of showing the distance of A strain to all the others
temp=sort.binary_dist_euclidean[sort.binary_dist_euclidean$strain1=="ECK0002-THRA",]
hist(temp$`Euclidean Distance`,xlab="Euclidean distance",
     main="Distance of ECK0002 to all the other ECKs")


