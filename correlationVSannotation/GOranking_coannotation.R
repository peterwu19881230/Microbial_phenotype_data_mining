#Use GO similarity as a distance metric and cumsum of co-annotated pairs on the y axis

str(strain1strain2_allAnnotations_allDistances)

anyCoAnnotation=ifelse(rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>0,
                       T,F)


#BP
##=====================================================================================================================
Wang_BP_anyCoAnnotation=data.frame(BP=strain1strain2_allAnnotations_allDistances$Wang_BP,CoAnnot=anyCoAnnotation)
Wang_BP_anyCoAnnotation=Wang_BP_anyCoAnnotation[!is.na(Wang_BP_anyCoAnnotation$BP),] #remove NA 
Wang_BP_anyCoAnnotation=Wang_BP_anyCoAnnotation[order(Wang_BP_anyCoAnnotation$BP,decreasing=T),]

#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr (BP) VS. Annot(any of 5)',
     xlab="high GO similarity (BP) -- rankings of strain pairs  -- low GO similarity (BP)",
     ylab="Cumsum of co-annotation",
     (1:dim(Wang_BP_anyCoAnnotation)[1])[samples],cumsum(Wang_BP_anyCoAnnotation$CoAnnot)[samples],
     type="l")
abline(a=0,b=sum(Wang_BP_anyCoAnnotation$CoAnnot)/dim(Wang_BP_anyCoAnnotation)[1],col="red",lty="dotted") #random expectation 
##=====================================================================================================================



#MF
##=====================================================================================================================
Wang_MF_anyCoAnnotation=data.frame(MF=strain1strain2_allAnnotations_allDistances$Wang_MF,CoAnnot=anyCoAnnotation)
Wang_MF_anyCoAnnotation=Wang_MF_anyCoAnnotation[!is.na(Wang_MF_anyCoAnnotation$MF),] #remove NA 
Wang_MF_anyCoAnnotation=Wang_MF_anyCoAnnotation[order(Wang_MF_anyCoAnnotation$MF,decreasing=T),]

#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr (MF) VS. Annot(any of 5)',
     xlab="high GO similarity (MF) -- rankings of strain pairs  -- low GO similarity (MF)",
     ylab="Cumsum of co-annotation",
     (1:dim(Wang_MF_anyCoAnnotation)[1])[samples],cumsum(Wang_MF_anyCoAnnotation$CoAnnot)[samples],
     type="l")
abline(a=0,b=sum(Wang_MF_anyCoAnnotation$CoAnnot)/dim(Wang_MF_anyCoAnnotation)[1],col="red",lty="dotted") #random expectation 
##=====================================================================================================================



#CC
##=====================================================================================================================
Wang_CC_anyCoAnnotation=data.frame(CC=strain1strain2_allAnnotations_allDistances$Wang_CC,CoAnnot=anyCoAnnotation)
Wang_CC_anyCoAnnotation=Wang_CC_anyCoAnnotation[!is.na(Wang_CC_anyCoAnnotation$CC),] #remove NA 
Wang_CC_anyCoAnnotation=Wang_CC_anyCoAnnotation[order(Wang_CC_anyCoAnnotation$CC,decreasing=T),]

#quick way to check the result
samples=1:4000 #note that this masks samples() function
plot(main='Corr (CC) VS. Annot(any of 5)',
     xlab="high GO similarity (CC) -- rankings of strain pairs  -- low GO similarity (CC)",
     ylab="Cumsum of co-annotation",
     (1:dim(Wang_CC_anyCoAnnotation)[1])[samples],cumsum(Wang_CC_anyCoAnnotation$CoAnnot)[samples],
     type="l")
abline(a=0,b=sum(Wang_CC_anyCoAnnotation$CoAnnot)/dim(Wang_CC_anyCoAnnotation)[1],col="red",lty="dotted") #random expectation 
##=====================================================================================================================