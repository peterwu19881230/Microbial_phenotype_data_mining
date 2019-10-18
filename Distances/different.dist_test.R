matrix=matrix(rnorm(100),ncol=10)

library(factoextra)
dist.mat.euclidean=get_dist(matrix,method="euclidean") ##get_dist is an enhanced version of dist()
dist.mat.maximum=get_dist(matrix,method="maximum")
dist.mat.manhattan=get_dist(matrix,method="manhattan")
dist.mat.canberra=get_dist(matrix,method="canberra")
dist.mat.binary=get_dist(matrix,method="binary")
dist.mat.minkowski=get_dist(matrix,method="minkowski")
dist.mat.spearman=get_dist(matrix,method="spearman")
dist.mat.kendall=get_dist(matrix,method="kendall")



require(gridExtra) ##This is required for creating side-by-side ggplot figures
p1=fviz_dist(dist.mat.euclidean) ##gives a ggplot figure
p2=fviz_dist(dist.mat.maximum)
p3=fviz_dist(dist.mat.manhattan)
p4=fviz_dist(dist.mat.canberra)
p5=fviz_dist(dist.mat.binary) ##This shows nothing because: https://stackoverflow.com/questions/23686028/how-the-command-distx-method-binary-calculates-the-distance-matrix
p6=fviz_dist(dist.mat.minkowski)
p7=fviz_dist(dist.mat.spearman)
p8=fviz_dist(dist.mat.kendall)

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=4)

