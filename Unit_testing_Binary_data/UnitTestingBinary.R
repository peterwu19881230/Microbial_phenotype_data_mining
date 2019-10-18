#Goal 1: find whether the clustering in R is reproducible by repeating the clustering using the same data:
binary_subset<-Binary_Data[1:100,1:100]
hclust<-vector("list",100) #c(,100) won't work. Also, "list" here is not a random text, it is a specified type
#Someone suggests to pre-allocate the size of the vector
# https://stackoverflow.com/questions/7815969/r-put-multiple-randomforest-objects-into-a-vector
for(i in 1:1000){
  dist<-dist(binary_subset,method="euclidean")
  hclust[[i]]<-hclust(dist,method="complete") #No idea what [[]] stands for
}

identical_clustering<-c()
for(i in 2:1000){
  identical_clustering<-c(identical_clustering,identical(hclust[[1]],hclust[[i]]))
}
sum(identical_clustering[TRUE]) #Should return 999 if they are all TRUE
#Result: they are all identical


