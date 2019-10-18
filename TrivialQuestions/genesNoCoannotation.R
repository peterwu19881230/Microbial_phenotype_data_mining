#How many gene pairs are not co-annotated in any of the annotation set above a particular similarity cutoff?

df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
TF=( rowSums(df[,-6])>=1 )
sum(TF)

cutoff=0.67

##no. of strains that are above the similarity cutoff:
sum((1-strain1strain2_allAnnotations_allDistances$pcc)>=0.67) #1962


##no. of strains that are co-annotated or not co-annotated

###co-annotated
sum((1-strain1strain2_allAnnotations_allDistances$pcc)>=0.67 & ( rowSums(df[,-6])>=1 )) #455


###not co-annotated
sum((1-strain1strain2_allAnnotations_allDistances$pcc)>=0.67 & ( rowSums(df[,-6])==0 )) #1507



1507/1962*100 #76.8% not co-annotated
