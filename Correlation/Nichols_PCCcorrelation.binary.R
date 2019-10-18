#Goal: Calculate all euclidean distances using binary data and create useful objects
##objects saved: binary_cor_strains, binary_sort.pcc.sql.NoIdent

#(!) Euclidean distance doesn't capture what pcc can capture for anti-correlation, so did some other distance method
#(!) Why is pcc not a good distance method for binary data?

#pcc for binary data
binary_cor_strains=cor(t(Binary_Data),use="pairwise.complete.obs",method="pearson") 
##According to the warning message: For sd=0, pcc would be NA 
## (Situations: 1. 2 rows only contain 0 and NA  2. Although 1 row has 1 or -1, but the position(s) for the other row is NA => ommited by use="pairwise.complete.obs")


## Euclidean distance for binary data
## The object defined in distances_and_hclusts.R:  binary_dist_euclidean


## Hamming distance (defined by sum of No. of unequal value for 2 rows)
## The object defined in distances_and_hclusts.R: 


##My opinion is that because there are very few missing values, pairwise.complete.obs won't cause significant difference

#Use No. to name the matrix. I can later retrieve ECK or other IDs from ECK_1st_table
rownames(binary_cor_strains)=1:3979
colnames(binary_cor_strains)=1:3979
save(binary_cor_strains,file="Data/binary_cor_strains.RData")




#Reorder binary_cor_strains into a SQL format (There should be a more elegant way: without converting distance obj to matrix)
##Sort the output of a dist obj: https://stackoverflow.com/questions/31591546/sorting-the-output-of-dist
m1=binary_cor_strains
obj_name="binary_sort.pcc.sql.NoIdent"

##--------------------
m1[upper.tri(m1)]=NA
diag(m1)=NA
##--------------------
##Can use these 2 lines if nothing=NA in the distance matrix (or correlation matrix). Although there are NA, it's ok since I want to remove them

library(reshape2)
m2=melt(m1) 
m2=m2[!is.na(m2$value),] ##Kick out all the NA
m2=m2[,c(2,1,3)] ##reorder the columns
names(m2)=c("strain1","strain2","Pearson.Correlation.Coefficient") ##name the columns
m2=m2[order(m2[,3],decreasing = T),] ##reorder by pcc
m2=m2[m2$strain1!=m2$strain2,]

binary_sort.pcc.sql.NoIdent=m2
save(binary_sort.pcc.sql.NoIdent,file=paste(obj_name,".RData",sep=""))
##Note: No. of id of strain1 and strain2 are both 3978 instead of 3979. This is correct 
##=> Easy to understand if we look at combination of: starin1=c(1,2,3,4) strain2=c(1,2,3,4). After taking out strain1==strain2 there are only 3 unique ids in strain1 and strain2

#Some stats
sd(binary_sort.pcc.sql.NoIdent$`Pearson.Correlation.Coefficient`)
hist(binary_sort.pcc.sql.NoIdent$`Pearson.Correlation.Coefficient`)



