#Goal: Do positive control using sig.3.pwy for the bottom up tree traversing experiment

#In order to pick up 3 genes from looking at the structure of the tree, I have to start with a smaller hclust
##Use objs in 3pathway.experiment.R (Copied so this file can be run directly)

library(dplyr)

dta=ECK.Pathway[order(ECK.Pathway$Data),]
sub.dta=dta[1:17,]

sig.3.pwy=All_Data[as.numeric(sub.dta$ids),]

##Use this to change labels to pathways
rownames(sig.3.pwy)=c("pw1-1","pw1-2","pw1-3","pw1-4","pw1-5","pw1-6","pw1-7","pw1-8","pw1-9","pw1-10",
                      "pw2-1","pw2-2",
                      "pw3-1","pw3-2","pw3-3","pw3-4","pw3-5")

dis=sig.3.pwy %>% pcc_dist ##pcc.dist is a function defined by me
hc=dis %>% hclust 
plot(hc)


##Nearest 3 genes
orderInTree=rbind(1:length(hc$labels),(1:length(hc$labels))[hc$order]) #1st row is the indices for original label, 2nd is their order in the tree
fractionBottomUp(hc,orderInTree[1,orderInTree[2,][3:5]]) ##3:5 mean the positions in the tree
noOfElementsBottomUp(hc,orderInTree[1,orderInTree[2,][3:5]])
fractionBottomUp(hc,orderInTree[1,orderInTree[2,][4:6]]) 
noOfElementsBottomUp(hc,orderInTree[1,orderInTree[2,][4:6]]) 

##Nearest 6 genes
fractionBottomUp(hc,orderInTree[1,orderInTree[2,][1:6]])
noOfElementsBottomUp(hc,orderInTree[1,orderInTree[2,][1:6]]) 
fractionBottomUp(hc,orderInTree[1,orderInTree[2,][12:17]])
noOfElementsBottomUp(hc,orderInTree[1,orderInTree[2,][12:17]])

##Nearest 11 genes
fractionBottomUp(hc,orderInTree[1,orderInTree[2,][7:17]]) 
noOfElementsBottomUp(hc,orderInTree[1,orderInTree[2,][7:17]]) 








