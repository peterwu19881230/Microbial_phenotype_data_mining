#Goal: Test the bottom up experiment. The algorithms are written at the top of functions defined in functions.R


library(dplyr)

#Test using a small portion of the data from ECK.Pathway:
dta=ECK.Pathway[order(ECK.Pathway$Data),]
sub.dta=dta[1:17,]

sig.3.pwy=All_Data[as.numeric(sub.dta$ids),]

##Use this to change labels to pathways
rownames(sig.3.pwy)=c("pw1-1","pw1-2","pw1-3","pw1-4","pw1-5","pw1-6","pw1-7","pw1-8","pw1-9","pw1-10",
                      "pw2-1","pw2-2",
                      "pw3-1","pw3-2","pw3-3","pw3-4","pw3-5")

dis=sig.3.pwy %>% pcc_dist ##pcc.dist is a function defined by me
hc=dis %>% hclust 

##Get the indices for the labels from str(hc)
rbind(hc$labels,hc$order)

##1st pathway
indices=1:10 ##This is the indices for the order of the labels in the original dataset
fractionBottomUp(hc,1:10)
noOfElementsBottomUp(hc,1:10)

##2nd pathway
indices=11:12
fractionBottomUp(hc,11:12)
noOfElementsBottomUp(hc,11:12)

##3rd pathway
indices=13:17
fractionBottomUp(hc,13:17)
noOfElementsBottomUp(hc,13:17)

##A more general way to do the whole experiment
lab=cbind(rownames(sig.3.pwy),c(rep("I",10),rep("II",2),rep("III",5))) %>% as.data.frame
result=bottomUpExp(hc,lab)

##Do its negative control
'
start.time=Sys.time()


j=c(2,5,10) # For No. of elements in these 3 pathways

##This block takes about 10min. It generates the fractions and so on. Resluts are in separate tables
count=0
avg.random_tables=list()
for(i in j){
  count=count+1
  random_tables=list()
  for (k in 1:10000){  ##How large should k be?
    indices=sample(rownames(sig.3.pwy),i)
    dat=data.frame(id=indices,tag=rep("random",i)) ##Give an arbitrary tag for the randomly picked up nodes
    random_tables[[k]]=bottomUpExp(hc,dat)
  }
  
  avg.fraction=c()
  avg.elements=c()
  NoOfElements=c()
  for(l in 1:length(random_tables)){
    avg.fraction[l]=random_tables[[l]]$"Avg Fraction"
    avg.elements[l]=random_tables[[l]]$"Average Elements in the Cluster"
    NoOfElements[l]=random_tables[[l]]$"No of Elements for the tag"
  }
  
  avg.avg.fraction=mean(avg.fraction)
  avg.avg.elements=mean(avg.elements)
  avg.NoOfElements=mean(NoOfElements)
  avg.random_tables[[count]]=data.frame(tag="random",avg.avg.fraction,avg.avg.elements,avg.NoOfElements,avg.NoOfElements/avg.avg.elements)
}
end.time=Sys.time()
time.taken=end.time-start.time
time.taken
test.avg.random_tables=avg.random_tables #assign to a new name so it does not mix with the one I do with the complete pathway experiment
#Object conversion from list to dataframe
test.avg.random_tables=rbind(test.avg.random_tables[[1]],test.avg.random_tables[[2]],test.avg.random_tables[[3]])
save(test.avg.random_tables,file="Data/test.avg.random_tables.RData")
'

 
### Create a vecotor of corresponding negative control (fraction obtained by radom sampling). Bind to auto.pwyBottomUpExp
randomFraction=sapply(result$`No of Elements for the tag`,FUN=function(number){
  index=which(number==test.avg.random_tables$avg.NoOfElements)
  return(test.avg.random_tables$avg.avg.fraction[index])
})

result$randomFraction=randomFraction

#Rename the experiment
three_pwy.bottomUp=result
save(three_pwy.bottomUp,file="Data/three_pwy.bottomUp.RData")

## Conclusion: the bottom-up experiment tells that the 3 pathways I chose as a toy example don't have a more clustered structure for each pathway
## Possible reasons: 1. bottomUp experiment has holes 2. Genes in each pathway of these 3 pathways are not more clutered


#----Note that the following doesn't contain negative controls----#

#Look at some specific pathways using a more manual way
##pick up the pathway: 1CMET2-PWY. Indices: c(800,1310,1461,1597,1628,1720,1809,2823,3422,3664)
indices=c(800,1310,1461,1597,1628,1720,1809,2823,3422,3664)
fractionBottomUp(hclust_pcc_complete,indices) #0.1 0.2 0.1 0.1 0.4 0.4 0.4 0.4 0.1 0.1
noOfElementsBottomUp(hclust_pcc_complete,indices) #Disappointing. It gives me 3979


##3 genes (Not all) that are in the pathway: GLYCLEAV-PWY (Nichols. fig 5B)
##ECK2898-GCVP, ECK2899-GCVH, ECK2900-GCVT: indices: c(1628,1720,1809)
indices=c(1628,1720,1809)
fractionBottomUp(hclust_pcc_complete,indices) #Perfect. All give me 1
noOfElementsBottomUp(hclust_pcc_complete,indices) #Perfect. All give me 3

##All genes from ALL-CHORISMATE-PWY
indices=(ECK.Pathway_table[ECK.Pathway_table$Pathway=="ALL-CHORISMATE-PWY",])$ids %>% as.numeric  
fractionBottomUp(hclust_pcc_complete,indices) #Gives the average = 0.07204861
noOfElementsBottomUp(hclust_pcc_complete,indices) #Disappointing. It gives me 3979


##All pathways (use s1.pwy and ECK.Pathway_table). 

##table for running the analysis (297 pathways. Note that from Pathways.col there are 299 pathways used)
temp=ECK.Pathway_table[ECK.Pathway_table$Pathway %in% s1.pwy,]
attach(temp)
unique.py=temp$Pathway %>% unique

all.py.fraction=list()
for(i in 1:length(unique.py)){
  indices=temp[Pathway==unique.py[i],]$ids %>% as.numeric
  print(i)
  print(indices)
  all.py.fraction[[i]]=fractionBottomUp(hclust_pcc_complete,indices)
}
names(all.py.fraction)=unique.py
##Result is in all.py.fraction. If the number is NaN (RHAMCAT-PWY ) that possibly means it's a 1 gene pathway (fraction=No. of found/total-1)


##Calculate the average:
avg.py.fraction=sapply(all.py.fraction,mean)  
avg.py.fraction.table=data.frame(avg.py.fraction)


##The following block takes 3min
all.py.NoOfElments=list()
for(i in 1:length(unique.py)){
  indices=temp[Pathway==unique.py[i],]$ids %>% as.numeric
  print(i)
  print(indices)
  all.py.NoOfElments[[i]]=noOfElementsBottomUp(hclust_pcc_complete,indices)
}
names(all.py.NoOfElments)=unique.py
##Result is in all.py.NoOfElments

##Calculate the average:
avg.py.NoOfElments=sapply(all.py.NoOfElments,mean)
avg.py.NoOfElments.table=data.frame(avg.py.NoOfElments)


##Bind a column that has the number of genes in the pathways
all.py.NoOfGenes=c()
for(i in 1:length(unique.py)){
  indices=temp[Pathway==unique.py[i],]$ids %>% as.numeric
  all.py.NoOfGenes[i]=length(indices)
}
names(all.py.NoOfGenes)=unique.py
##Result is in all.py.NoOfGenes

##Create a summarizing table
pwy.bottomToTop.table=cbind(avg.py.fraction.table,avg.py.NoOfElments.table,all.py.NoOfGenes,all.py.NoOfGenes/avg.py.NoOfElments.table)
colnames(pwy.bottomToTop.table)=c("fraction","NoOfElements","NoOfGenes","NoOfGenes/NoOfElements")
###For those who have NoOfGenes=1 and fraction=NA and so on: Only 1 gene is found among all ECKs by using the 299 pathway names retrieved from Pathways.col 
###=> Can be verified by looking at EcoCycID.Pwys



mean(pwy.bottomToTop.table$fraction[pwy.bottomToTop.table$fraction!=999]) ##The mean of the fraction is 0.2646306
detach(temp)

