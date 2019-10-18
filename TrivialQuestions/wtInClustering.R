# Goal: determine how scattered the wildtypes are using quantitative scores to cluster 
# (Wildtypes are defined by the new cutoff=-3.463 illustrated in the file phenptypeCutoff.R) 

## Get the wildtype label (row indices) from All_Data
new.binary=All_Data
new.binary=ifelse(new.binary<=-3.463 | new.binary>=3.463,1,0)

TFvector=apply(new.binary,1,FUN=function(x){
  row.sum=sum(x,na.rm=T)
  ifelse(row.sum==0,T,F)
})

indices=(1:3979)[TFvector]

## Draw the dendrogram from the tree and color the wiltypes (Haven't found a way to remove the labels)
library(dendextend); library(dplyr)
hclust_pcc_complete %>% as.dendrogram %>% set("by_labels_branches_col", value = indices) %>% plot

rm(TFvector,indices) ## remove the unnecessary objects