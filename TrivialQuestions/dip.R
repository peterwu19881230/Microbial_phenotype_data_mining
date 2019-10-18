#Goal: Figure out what the dip is in CorrelationVSAnnotation experiment. I use the experiment that combines pwy and ptcomplex

if (!exists("pccPwyPcom")) {
  
  # Quantitative data preparation pwy
  start.time = Sys.time()
  strain1strain2.pcc.pwyTF = merge(abs.sort.pcc.sql.NoIdent, strain1strain2.samePWY_Annot, by = c("strain1", "strain2"))
  end.time = Sys.time()
  end.time - start.time  #Time difference of 2.562144 mins
  
  ## pcomplex take abs of pcc
  strain1strain2.pcc.pcomTF = strain1strain2.pcc.samePTcomplex_Annot
  strain1strain2.pcc.pcomTF$Pearson.Correlation.Coefficient = abs(strain1strain2.pcc.pcomTF$Pearson.Correlation.Coefficient)
  
  
  ## pwy+pcomplex
  start.time = Sys.time()
  strain1strain2.pcc.samePwyAnnot.samePcomAnnot = merge(strain1strain2.pcc.pwyTF, 
                                                        strain1strain2.pcc.pcomTF, by = c("strain1", "strain2"))
  end.time = Sys.time()
  end.time - start.time  #Time difference of 1.128728 mins
  
  
  pcc_TF = apply(strain1strain2.pcc.samePwyAnnot.samePcomAnnot[, c("bi.pwy.annot", 
                                                                   "bi.ptcomplex.annot")], 1, FUN = function(row) {
                                                                     if (row[1] == 0 & row[2] == 0) {
                                                                       return(0)
                                                                     } else {
                                                                       1
                                                                     }
                                                                   })
  
  ## A table that is: strain1 - strain2 - pcc - PwyPcom_TF
  pccPwyPcom = cbind(strain1strain2.pcc.samePwyAnnot.samePcomAnnot[, 1:3], 
                     pcc_TF)
  pccPwyPcom = pccPwyPcom[order(pccPwyPcom$Pearson.Correlation.Coefficient.x, 
                                decreasing = T), ]  ##sort by pcc
  save(pccPwyPcom, file = "Data/pccPwyPcom.RData")
}

if (!exists("hammingPwyPcom")) {
  # Binary data preparation pwy
  start.time = Sys.time()
  strain1strain2.hammingNo0.samePWY_Annot = merge(binary_coef_hamming_No0,strain1strain2.samePWY_Annot, by = c("strain1", "strain2"))
  end.time = Sys.time()
  end.time - start.time  #Time difference of 1.25747 mins
  
  ## pcomplex
  start.time = Sys.time()
  strain1strain2.hammingNo0.samePTcomplex_Annot = merge(binary_coef_hamming_No0, strain1strain2.pcc.samePTcomplex_Annot[, -3], by = c("strain1", "strain2"))
  end.time = Sys.time()
  end.time - start.time  #Time difference of 1.246467 mins
  
  ## pwy+pcomplex
  start.time = Sys.time()
  strain1strain2.hammingNo0.samePwyAnnot.samePcomAnnot = merge(strain1strain2.hammingNo0.samePWY_Annot, strain1strain2.hammingNo0.samePTcomplex_Annot, by = c("strain1", "strain2"))
  end.time = Sys.time()
  end.time - start.time  #Time difference of 42.60484 secs
  
  hammingNo0_TF = apply(strain1strain2.hammingNo0.samePwyAnnot.samePcomAnnot[, 
                                                                             c("bi.pwy.annot", "bi.ptcomplex.annot")], 1, FUN = function(row) {
                                                                               if (row[1] == 0 & row[2] == 0) {
                                                                                 return(0)
                                                                               } else {
                                                                                 1
                                                                               }
                                                                             })
  
  
  ## A table that is: strain1 - strain2 - hamNo0 - PwyPcom_TF
  hammingPwyPcom = cbind(strain1strain2.hammingNo0.samePwyAnnot.samePcomAnnot[,1:3], hammingNo0_TF)
  hammingPwyPcom = hammingPwyPcom[order(hammingPwyPcom$`Hamming Distance.x`, 
                                        decreasing = T), ]  ##sort by hamNo0
  save(hammingPwyPcom, file = "Data/hammingPwyPcom.RData")
}


# The following fig looks too much the same as pwy fig alone. Below is my
# verification that Fig.1 is correct:
# sum(strain1strain2.pcc.samePwyAnnot.samePcomAnnot$bi.pwy.annot) #7787
# sum(strain1strain2.pcc.samePwyAnnot.samePcomAnnot$bi.ptcomplex.annot)
# #1221 sum(pcc_TF) #8782

## sum(strain1strain2.hammingNo0.samePwyAnnot.samePcomAnnot$bi.pwy.annot)
## #7787
## sum(strain1strain2.hammingNo0.samePwyAnnot.samePcomAnnot$bi.ptcomplex.annot)
## #1221 sum(hammingNo0_TF) #8782



TF = pccPwyPcom$pcc_TF
TF2 = hammingPwyPcom$hammingNo0_TF
fraction = cumsum(TF)/(1:length(TF))
fraction2 = cumsum(TF2)/(1:length(TF2))
coef=pccPwyPcom$Pearson.Correlation.Coefficient.x
coef2=hammingPwyPcom$`Hamming Distance.x`


plot(main = "Fig. 1(b): First 10000 coef V.S. pwy or pcom annotation in Nichols'\n (cummulative for coef)", 
     1:10000, fraction[1:10000], type = "l", xlab = "Order of the decreasing coef", 
     ylab = "No. of the same pwy or pcom annotation / No. of coef until the  cutoff coef", 
     ylim = c(0, 1))
text(c(1, 2000, 4000, 6000, 8000, 10000), fraction[c(1, 2000, 4000, 6000, 8000, 
                                                     10000)], pos = 3, labels = round(c(coef[1], coef[2000], coef[4000], coef[6000], 
                                                                                        coef[8000], coef[10000]), 4), cex = 0.5)  #'3' in pos=3 means 'above' 
lines(1:10000, fraction2[1:10000], type = "l", col = "red", ylim = c(0, 1))
text(c(1, 2000, 4000, 6000, 8000, 10000), fraction2[c(1, 2000, 4000, 6000, 8000, 
                                                      10000)], pos = 3, labels = round(c(coef2[1], coef2[2000], coef2[4000], coef2[6000], 
                                                                                         coef2[8000], coef2[10000]), 4), cex = 0.5, col = "red")  #'3' in pos=3 means 'above'

#I made the 2 lines by eye-balling the dip
abline(v=30,col="blue")
abline(v=300,col="blue") 


plot(main = "Fig. 1(c): First 1000 coef V.S. pwy or pcom annotation in Nichols' \n(cummulative for coef)", 
     1:1000, fraction[1:1000], type = "l", xlab = "Order of the decreasing coef", 
     ylab = "No. of the same pwy or pcom annotation / No. of coef until the  cutoff coef", 
     ylim = c(0, 1))
text(c(1, 200, 400, 600, 800, 1000), fraction[c(1, 200, 400, 600, 800, 1000)], pos = 3, 
     labels = round(c(coef[1], coef[200], coef[400], coef[600], coef[800], coef[1000]), 
                    4), cex = 0.5)  #'3' in pos=3 means 'above'
lines(1:1000, fraction2[1:1000], type = "l", col = "red", ylim = c(0, 1))
text(c(1, 200, 400, 600, 800, 1000), fraction2[c(1, 200, 400, 600, 800, 1000)], pos = 3, 
     labels = round(c(coef2[1], coef2[200], coef2[400], coef2[600], coef2[800], coef2[1000]), 
                    4), cex = 0.5, col = "red")  #'3' in pos=3 means 'above'

# I made the 2 lines by eye-balling the dip
abline(v=30,col="blue")
abline(v=300,col="blue") 

dip_table_pcc=pccPwyPcom[30:300,]
dip_strains=unique(c(unique(dip_table_pcc$strain1),unique(dip_table_pcc$strain2)))
dip_strain_table=ECK_1st_table[dip_strains,]
##write.table(dip_table_pcc,file="dip_table_pcc.txt",row.names = F,quote=F)
##write.table(dip_strain_table,file="dip_strain_table.txt",row.names = F,quote=F)

# Add columns for names and reorder the columns
dip_table_pcc_2=dip_table_pcc
dip_table_pcc_2$strain1_gene_names=ECK_1st_table$associated_gene_names[dip_table_pcc$strain1]
dip_table_pcc_2$strain2_gene_names=ECK_1st_table$associated_gene_names[dip_table_pcc$strain2]
dip_table_pcc_2=dip_table_pcc_2[,c("strain1","strain1_gene_names","strain2","strain2_gene_names",
                               "Pearson.Correlation.Coefficient.x","pcc_TF"
                               )]

##write.table(dip_table_pcc_2,file="dip_table_pcc_2.txt",row.names = F,quote=F) 
##Note: This file still contains genes involved in the same pwy or ptcomplex

# Remove those genes that are involved in the same pwy or ptcomplex
dip_table_pcc_3=dip_table_pcc_2[,c(2,4,5,6)] ##remove the strain ids (I don't need them here)
dip_table_pcc_3=dip_table_pcc_3[dip_table_pcc_3$pcc_TF==0,] ## Remove genes that are involved in the same pwy or ptcomplex
dip_table_pcc_3=dip_table_pcc_3[,-4] #I don't need the last column now
##write.table(dip_table_pcc_3,file="dip_table_pcc_3.txt",row.names = F,quote=F)

#Note: the files generated above are now in the Data directory
