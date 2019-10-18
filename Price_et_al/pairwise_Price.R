#construct: strain1 - strain2 (strains in Price's that are also in Nichols') - pcc - pwy co-annotation - pcomplex co-annotation - ....etc

load("Data/Nich_Price_quantitative.RData")
Price_quantitative=Nich_Price_quantitative[,325:486]
Price_pairwisePCC_longTable=( 1-abs(cor(t(Price_quantitative))) ) %>% as.dist %>% meltANDsort_dist #the 3rd column is the pcc based distance
names(Price_pairwisePCC_longTable)=c("strain1","strain2","1-|PCC|")
Price_pairwisePCC_longTable$strain1=as.numeric(Price_pairwisePCC_longTable$strain1)
Price_pairwisePCC_longTable$strain2=as.numeric(Price_pairwisePCC_longTable$strain2)


Price_pairwisePCC_longTable_coannotations=dplyr::left_join(Price_pairwisePCC_longTable,
                       strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","Pwy","pcomplex","operon","regulator","kegg_modules")],
                       by=c("strain1","strain2"))

###anyNA(Price_pairwisePCC_longTable_coannotations$Pwy) #wrong order of strain1-strain2 in the Price table generates NAs

##swap and rejoin (using pwy coannotation column alone should be enough):
Price_pairwisePCC_longTable_coannotations[is.na(Price_pairwisePCC_longTable_coannotations$Pwy),
                                          c("strain1","strain2")]=
Price_pairwisePCC_longTable_coannotations[is.na(Price_pairwisePCC_longTable_coannotations$Pwy),
                                            c("strain2","strain1")]

Price_pairwisePCC_longTable_coannotations=dplyr::left_join(Price_pairwisePCC_longTable_coannotations[,c("strain1","strain2","1-|PCC|")],
                                                           strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","Pwy","pcomplex","operon","regulator","kegg_modules")],
                                                           by=c("strain1","strain2"))

##all the NAs are gone
###anyNA(Price_pairwisePCC_longTable_coannotations$Pwy)
###anyNA(Price_pairwisePCC_longTable_coannotations$pcomplex)
###anyNA(Price_pairwisePCC_longTable_coannotations$operon)
###anyNA(Price_pairwisePCC_longTable_coannotations$regulator)
###anyNA(Price_pairwisePCC_longTable_coannotations$kegg_modules)
  
#save(Price_pairwisePCC_longTable_coannotations,file="Data/Price_pairwisePCC_longTable_coannotations.RData")



TF=( Price_pairwisePCC_longTable_coannotations[,c("Pwy","pcomplex","operon","regulator","kegg_modules")] %>% rowSums ) ==5
Price_pairwisePCC_cumsum=Price_pairwisePCC_longTable_coannotations[,c("strain1","strain2","1-|PCC|")] %>% cbind(TF)
Price_pairwisePCC_cumsum$cumsum_all=cumsum(Price_pairwisePCC_cumsum$TF)
Price_pairwisePCC_cumsum=Price_pairwisePCC_cumsum[,-4] #remove the TF column

#save(Price_pairwisePCC_cumsum,file="Data/Price_pairwisePCC_cumsum.RData")


