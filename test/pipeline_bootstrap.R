#pipeline to do bootstrap test
##Using this ref is too slow so I don't consider it:https://stats.stackexchange.com/questions/20701/computing-p-value-using-bootstrap-with-r (Note: some code in STAT646 by Moumita also has it.)


##data prep. Code reused from: AvgPCC_pwy$ptcomplex.R 
##===================================================
distance_column="pcc"

##Pwy
coAnnotated= ( strain1strain2_allAnnotations_allDistances$Pwy==1 )
pwy_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated] #I am converting the |PCC| based distance back to just |PCC|

##pcomplex
coAnnotated= ( strain1strain2_allAnnotations_allDistances$pcomplex==1 )
pcomplex_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]


##Pwy and pcomplex 
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex")])==2 )
pwyANDpcomplex_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]


#all_annotSet
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5 )
all_annotSet_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]



# all
All=1-strain1strain2_allAnnotations_allDistances[[distance_column]]
##===================================================

#(!) (The following isn't complete: Neither the code nor the concept)

#Ref for the concept: https://campus.datacamp.com/courses/statistical-thinking-in-python-part-2/introduction-to-hypothesis-testing?ex=10

#determine bootstrap pvalues

#dat1=All
#dat2=pwy_abs_pcc
dat1=c(1,1,1,1,1.1,1.1,1.2)
dat2=c(5,6,4.1,5,5,5)

n_of_boot=5000 #no. of boostrap sampling

dat1_all=matrix(rep(dat1,n_of_boot),ncol=n_of_boot,byrow = F)
dat2_all=matrix(rep(dat2,n_of_boot),ncol=n_of_boot,byrow = F)

resample_mean=function(values){
  mean(sample(values,replace=T))
}


#compute bootstrap p-val
mean( 
  ( apply(dat2_all,2,FUN=resample_mean) - apply(dat1_all,2,FUN=resample_mean) ) <= 0
) #0















