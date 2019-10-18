#Goal:Fig S1 (pathway heatmap)

##Get the pathway Unique ID that Nichols used (retrieve from pathwayscol instead of from Nichols to get the correct No.: 299)
uniq=unique(EcoCycID.Pwys$`Pathway ID`)

s1.pwy=c()
j=1
for(i in 1:length(uniq)){
  indices=which(EcoCycID.Pwys$`Pathway ID` %in% uniq[i]) 
  if(length(indices)>2){
    s1.pwy[j]=EcoCycID.Pwys$`Pathway ID`[indices[1]] ##only need 1 index for retrieving the pathway ID
    j=j+1
  } 
} ##this gives 300 pathways

##drop the tRNA charging pwy which contains 109 genes
s1.pwy=s1.pwy[s1.pwy!="TRNA-CHARGING-PWY"] ##then s1.pwy contains exactly 299 pwys
save(s1.pwy,file="Data/s1.pwy.RData")


##Construct the pathway signature matrix
library(dplyr)
pwy.sig=data.frame()
for(i in 1:length(s1.pwy)){
    indices=which(ECK.Pathway$Data %in% s1.pwy[i]) ### the indices for grabbing the ECK ids from ECK.Pathway (ids were assigned according to rownames of allData.csv )
    ids=ECK.Pathway$ids[indices] %>% as.numeric ###ECK ids for a pathway
    if(identical(ids,integer(0))==F) pwy.sig=rbind(pwy.sig,colMeans(All_Data[ids,],na.rm=T))
        
}
colnames(pwy.sig)=colnames(All_Data)
rownames(pwy.sig)=s1.pwy
pwy.sig=as.matrix(pwy.sig) ##I don't know how to initiate pwy.sig using matrix so I started with data frame and converted it back here

##From below I see I have to remove some NAN values
##Heatmap(pwy.sig,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,col=colorRamp2(c(-10,0,10),c("black","black","black")))
pwy.sig=pwy.sig[-c(253,264),] ##I am just removing the 2 rows which contains only NAN so that Heatmap() can run at least create the correlation matrix 


dist.pwy.sig=cor(t(pwy.sig),use="pairwise.complete.obs", method="pearson")

NAs=which(is.na(dist.pwy.sig),arr.ind=T) ## Where the cor matrix contains NAs
###How should I deal with these NAs?

library(ComplexHeatmap)
library(circlize)
Heatmap(dist.pwy.sig,show_row_names = F,show_column_names = F,col=colorRamp2(c(-1,-0.67,-0.33,0,0.33,0.67,1),c("cyan","#1785A2","#144F57","black","#6E701A","#A7AA29","yellow"))) ##When running a subset it is ok. I think there are some NA values
##The default clustering method is complete linkage
##The reuslt looks kind of similar to fig.S1 but there are lots of differences

##The color code was done by eye-picking from: http://htmlcolorcodes.com  v.s. original figS1. It seems the colors are continuous, not discrete (because it looks the same as the one that uses only (c(-1,0,1),c("cyan","black","yellow")))



## 2 things that I am not sure about the result: 
##1. Did they use na.rm for colMeans of the pathways?
##2. Did they use: pairwise.complete.obs for cor()?





