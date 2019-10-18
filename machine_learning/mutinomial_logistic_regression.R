#Ref:https://statisticsbyjim.com/regression/choosing-regression-analysis/
#Ref:https://en.wikipedia.org/wiki/Multinomial_logistic_regression


load("Data/annot_20_clusters_subset_hclust_hamming.RData")

annot_20_clusters_subset

dat=cbind(annot_20_clusters_subset,All_Data_NAimputed[names(annot_20_clusters_subset),])


## I don't fully understand multinomial-logistic-regression, the following is just a test based on:
## Ref:https://stats.idre.ucla.edu/r/dae/multinomial-logistic-regression/

out=multinom(annot_20_clusters_subset ~ dat[,1] + dat[,2] + dat[,3])

predict(out,newdata=data.frame(1,2,3),"probs")
### Error in Y[keep, ] <- Y1 : NAs are not allowed in subscripted assignments
### In addition: Warning message:
### 'newdata' had 1 row but variables found have 147 rows 
