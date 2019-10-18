#Goal: Calculate all correlation coefficients using all Nichols' data and create useful objects
##objects saved: cor_strains, sort.pcc.sql.NoIdent, cor_conditions, sort.ConditionPCC.sql.NoIdent 


#A new strainCC object
'
cor_strains=cor(t(All_Data),use="pairwise.complete.obs",method="pearson") ##This is the question for Nasos (whether they used pairwise.complete.obs)

#Use No. to name the matrix. I can later retrieve ECK or other IDs from ECK_1st_table
rownames(cor_strains)=1:3979
colnames(cor_strains)=1:3979

save(cor_strains,file="Data/sourced/cor_strains.RData")
'
##My opinion is that because there are very few missing values, pairwise.complete.obs won't cause significant difference

### No. of NA
sum(is.na(as.numeric(as.matrix(All_Data)))) #2234. NAs/total=2234/1289196=0.001732863 ~ 0.2%





#Reorder cor_strains into a SQL format
##Sort the output of a dist obj: https://stackoverflow.com/questions/31591546/sorting-the-output-of-dist
m1=cor_strains
m1[upper.tri(m1)]=0 ## This can be done when no value in the pcc is 0
library(reshape2)
m2=subset(melt(m1), value!=0)
m2=m2[,c(2,1,3)] ##reorder the columns
names(m2)=c("strain1","strain2","Pearson.Correlation.Coefficient") ##name the columns
m2=m2[order(m2[,3],decreasing = T),] ##reorder by PCC 
sort.pcc.sql=m2; rm(m2)
sort.pcc.sql.NoIdent=sort.pcc.sql[sort.pcc.sql$strain1!=sort.pcc.sql$strain2,]
##save(sort.pcc.sql.NoIdent,file="Data/sourced/sort.pcc.sql.NoIdent.RData")
##Note: No. of id of strain1 and strain2 are both 3978 instead of 3979. This is correct 
##=> Easy to understand if we look at combination of: starin1=c(1,2,3,4) strain2=c(1,2,3,4). After taking out strain1==strain2 there are only 3 unique ids in strain1 and strain2

#Some stats
sd(cor_strains)
hist(cor_strains)

#From IRIS paper (2016): their comparison for the correlations (fig. 2D)
IRIS_Nichols.cor_strains=read.csv("Data/opacity.s.scores.cor.vs.size.s.scores.cor.csv")
OriginalCC=IRIS_Nichols.cor_strains$correlation.published
IntOpaCC=IRIS_Nichols.cor_strains$correlation.opacity

#Density scatter plot
smoothScatter(OriginalCC,IntOpaCC,colramp = colorRampPalette(c("white", blues9))) 
smoothScatter(OriginalCC,IntOpaCC,colramp = colorRampPalette(c("black", blues9))) 
##Looks weird because the center is white => 
Original.Int=cbind(OriginalCC,IntOpaCC)
between=Original.Int[(OriginalCC>-0.2&OriginalCC<0.2)&(IntOpaCC>-0.2&IntOpaCC<0.2)] #0 elements
between=Original.Int[(OriginalCC>=-0.2&OriginalCC<=0.2)&(IntOpaCC>=-0.2&IntOpaCC<=0.2)] # 14 elements

##Remove diagonal elements. This is just my guess about what they have done. I don't know which sd() they use so I just pick up sd(OriginalCC)
sd.Ori=sd(OriginalCC)
sd.Int=sd(IntOpaCC)
smoothScatter(OriginalCC[abs(OriginalCC-IntOpaCC)<=3*sd(OriginalCC)],IntOpaCC[abs(OriginalCC-IntOpaCC)<=3*sd(OriginalCC)],colramp = colorRampPalette(c("white", blues9)))

## Points above 3 sd
sum(abs(OriginalCC-IntOpaCC)>3*sd(OriginalCC)) ## 108 points that are above 3 sd
sum((OriginalCC-IntOpaCC)>3*sd(OriginalCC)) ##green: 94 
sum((IntOpaCC-OriginalCC)>3*sd(OriginalCC)) ##red 14
###It seems to me the paper is BSing because the numbers don't match theirs. They have way more red and green dots

