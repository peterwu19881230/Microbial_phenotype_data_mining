#Condition correlation coefficient object



cor_conditions=cor(All_Data,use="pairwise.complete.obs",method="pearson") ##This is the question for Nasos (whether they used pairwise.complete.obs)
#save(cor_conditions,file="Data/sourced/cor_conditions.RData")

#Reorder cor_conditions into a SQL format (There should be a more elegant way: without converting distance obj to matrix)
##Sort the output of a dist obj: https://stackoverflow.com/questions/31591546/sorting-the-output-of-dist

sort.ConditionPCC.sql.NoIdent=meltANDsort_dist(cor_conditions)
#save(sort.ConditionPCC.sql.NoIdent,file="Data/sourced/sort.ConditionPCC.sql.NoIdent.RData")
##Note: No. of id of strain1 and strain2 are both 3978 instead of 3979. This is correct 
##=> Easy to understand if we look at combination of: starin1=c(1,2,3,4) strain2=c(1,2,3,4). After taking out strain1==strain2 there are only 3 unique ids in strain1 and strain2

#Some stats
sd(cor_conditions)
hist(cor_conditions)
