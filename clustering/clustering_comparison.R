#Clustering of original data
mat.fig.4A.strain<-All_Data[which(ECK_1st_table$sorted_ECK_missing_gene_names_added %in% c("ECK0457-ACRA","ECK0456-ACRB","ECK1524-MARA","ECK1525-MARB","ECK1523-MARR")),]
dist.mat.fig.4A.strain<-dist(mat.fig.4A.strain,method="euclidean")
hclust.mat.fig.4A.strain<-hclust(dist.mat.fig.4A.strain,method="complete")

'

#hclust content:

structure(list(
merge = structure(c(-4L, -3L, -1L, 2L, -5L, 1L, -2L, 3L), .Dim = c(4L, 2L)), 
height = c(26.0838443996607, 46.1777521326489, 81.6352250717119, 173.772157800603), 
order = c(3L, 4L, 5L, 1L, 2L), 
labels = c("1278" "2432" "2693" "2775" "3965"), 
method = "complete", 
call = hclust(d = dist.mat.fig.4A.strain, method = "complete"), 
dist.method = "euclidean"), .Names = c("merge", "height", "order", "labels", "method", "call", "dist.method"), class = "hclust")

'


##plot(hclust_original,main = "Original (euclidean + complete)")



#Clustering of the binary data: 5 genotypes from fig.4A using all 324 conditions 
mat.binary.fig.4A.strain<-Binary_Data[which(ECK_1st_table$sorted_ECK_missing_gene_names_added %in% c("ECK0457-ACRA","ECK0456-ACRB","ECK1524-MARA","ECK1525-MARB","ECK1523-MARR")),]
dist.mat.binary.fig.4A.strain<-dist(mat.binary.fig.4A.strain,method="euclidean")
hclust.mat.binary.fig.4A.strain<-hclust(dist.mat.binary.fig.4A.strain,method="complete")
dendro.mat.binary.fig.4A.strain<-as.dendrogram(hclust.mat.binary.fig.4A.strain)


'
This part tries to add more genotypes and re-cluster

#Add more genotypes and re-do the clustering
mat.binary.fig.4A.strain<-Binary_Data[which(ECK_1st_table$sorted_ECK_missing_gene_names_added %in% c("ECK0457-ACRA","ECK0456-ACRB","ECK1524-MARA","ECK1525-MARB","ECK1523-MARR")),]
matrix<-mat.binary.fig.4A.strain

##Excluding the existing 5 genotypes from the binary data for latter use:
index=which(ECK_1st_table$sorted_ECK_missing_gene_names_added %in% c("ECK0457-ACRA","ECK0456-ACRB","ECK1524-MARA","ECK1525-MARB","ECK1523-MARR"))
restOfBinary=Binary_Data[-index,]
rm(index)

##Add 1 more genotype each time and store them in a list called hclust_binary_extended_dendro

hclust_binary_extended_dendro<-list()
for(i in 1:20){  ##Change this number to add more genotypes
matrix<-rbind(matrix,Binary_Data[1:i,]) 
dendro<-matrix %>% dist(,method="euclidean") %>% hclust(,method="complete") %>% as.dendrogram
hclust_binary_extended_dendro[[i]]<-dendro
}

'


'
This part compares the differences of trees after addition of more genotypes

##Compare the dendrogram (Use.edge.length to ignore the heights of the branches of the tree. Described in: https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html
##=>"Global Comparison of two (or more) dendrograms")

###Prune the trees so they only contain the original 5 genes .If I put the argument name in front of each argument, it complains about not finding them (blah blah is missing with no default). Very strage. 
DendroPruningBy5Genes<-lapply(hclust_binary_extended_dendro,intersect_trees,dend2=dendro5genes)

###Pairwise comparison to the original clustering of 5 genes. The result is a sequence of booleans
TreeComparisons<-sapply(DendroPruningBy5Genes,all.equal,dend2=dendro5genes,use.edge.length = FALSE)

#No. of identical trees and non-identical trees
sink("fig4A_tree_compare.txt")
print("No. of identical trees:")
sum(TreeComparisons[TreeComparisons==TRUE])
print("No. of non-identical trees:")
sum(TreeComparisons[TreeComparisons==FALSE])

'


'
#What the data look like:
##Original
H1<-Heatmap(name = "OriginalScores",mat.fig.4A.strain,cluster_rows=FALSE,cluster_columns=FALSE,show_column_names = FALSE)
##Binary
H2<-Heatmap(name = "BinaryScores",mat.binary.fig.4A.strain,cluster_rows=FALSE,cluster_columns=FALSE,show_column_names = FALSE)
H1+H2
'



#I had an error running this. Don't know why:
'
library("pvclust")
  
  #Original
> result1<-pvclust(t(mat.fig.4A.strain),method.dist = "euclidean",method.hclust = "complete",nboot=1000)
Bootstrap (r = 0.4)... Done.
Bootstrap (r = 0.6)... Done.
Bootstrap (r = 0.6)... Done.
Bootstrap (r = 0.8)... Done.
Bootstrap (r = 0.8)... Done.
Bootstrap (r = 1.0)... Done.
Bootstrap (r = 1.0)... Done.
Bootstrap (r = 1.2)... Done.
Bootstrap (r = 1.2)... Done.
Bootstrap (r = 1.4)... Done.
Error in solve.default(crossprod(X, X/vv)) : 
  Lapack routine dgesv: system is exactly singular: U[2,2] = 0
In addition: There were 11 warnings (use warnings() to see them)


#Binary
> result2<-pvclust(t(mat.binary.fig.4A.strain),method.dist = "euclidean",method.hclust = "complete",nboot=1000)
'
