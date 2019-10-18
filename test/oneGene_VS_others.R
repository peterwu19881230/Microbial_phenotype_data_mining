#Look at how a gene in Nichols' correlate with others. Will I be able to gain insight by these other 3978 genes' pcc?

temp=strain1strain2_allAnnotations_allDistances[strain1strain2_allAnnotations_allDistances$strain1==1,c("strain1","strain2","pcc")]
temp$pcc=1-temp$pcc

str(temp)
hist(temp$pcc)
summary(temp$pcc)