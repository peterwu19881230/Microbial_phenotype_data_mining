#Goal: See if different distances can help sort the relavence of ECKs

#Sort the output of a dist obj: https://stackoverflow.com/questions/31591546/sorting-the-output-of-dist


sort.dist.All_Data_NAimputed.euclidean=meltANDsort_dist(dist.All_Data_NAimputed.euclidean)
sort.dist.All_Data_NAimputed.maximum=meltANDsort_dist(dist.All_Data_NAimputed.maximum)
sort.dist.All_Data_NAimputed.manhattan=meltANDsort_dist(dist.All_Data_NAimputed.manhattan)
sort.dist.All_Data_NAimputed.canberra=meltANDsort_dist(dist.All_Data_NAimputed.canberra)
sort.dist.All_Data_NAimputed.binary=meltANDsort_dist(dist.All_Data_NAimputed.binary)
sort.dist.All_Data_NAimputed.minkowski=meltANDsort_dist(dist.All_Data_NAimputed.minkowski)
sort.dist.All_Data_NAimputed.spearman=meltANDsort_dist(dist.All_Data_NAimputed.spearman)
#(!)For kendall I wasn't even able to get the distance matrix (too slow to have finished on my Mac)

names(sort.dist.All_Data_NAimputed.euclidean)=c("strain1","strain2","euclidean")
names(sort.dist.All_Data_NAimputed.maximum)=c("strain1","strain2","maximum")
names(sort.dist.All_Data_NAimputed.manhattan)=c("strain1","strain2","manhattan")
names(sort.dist.All_Data_NAimputed.canberra)=c("strain1","strain2","canberra")
names(sort.dist.All_Data_NAimputed.binary)=c("strain1","strain2","binary")
names(sort.dist.All_Data_NAimputed.minkowski)=c("strain1","strain2","minkowski")
names(sort.dist.All_Data_NAimputed.spearman)=c("strain1","strain2","spearman")

save( #Saving these files take a while...
  sort.dist.All_Data_NAimputed.euclidean,
  sort.dist.All_Data_NAimputed.maximum,
  sort.dist.All_Data_NAimputed.manhattan,
  sort.dist.All_Data_NAimputed.canberra,
  sort.dist.All_Data_NAimputed.binary,
  sort.dist.All_Data_NAimputed.minkowski,
  sort.dist.All_Data_NAimputed.spearman,
  file="Data/meltedANDsorted_distances.RData"
)

#example of showing the distance of A strain to all the others (Euclidean)
temp=sort.dist.All_Data_NAimputed.euclidean[sort.dist.All_Data_NAimputed.euclidean$strain1==1,]
hist(temp$`Euclidean Distance`,xlab="Euclidean distance",
     main="Distance of 1 to all the other ECKs",breaks=100)



