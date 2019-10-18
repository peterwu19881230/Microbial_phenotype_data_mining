df=strain1strain2_allAnnotations_allDistances
pcc_all=1-df$pcc
pcc_pwy=1-df$pcc[df$Pwy==1] #sample 1
pcc_pcom=1-df$pcc[df$pcomplex==1] #sample 2
pcc_pwyANDpcom=1-df$pcc[df$Pwy==1 & df$pcomplex==1] #sample 3
pcc_5Annot=1-df$pcc[rowSums(df[,c("Pwy","pcomplex","regulator","operon","kegg_modules")])==5] #sample 4

#save(pcc_all,pcc_pwy,pcc_pcom,pcc_pwyANDpcom,pcc_5Annot,file="pcc.RData")

##pcc data are uploaded to my Github: "https://github.com/peterwu19881230/test_repo/blob/master/pcc.RData"

#load("pcc.RData") 

my_transform=function(x){ x/(1-x) } 

H0=my_transform(pcc_all) #I use all pcc to approximate H0. This doesn't follow normal
#(!) I cannot do |r1|/(1-|r1|) - |r2|/(1-|r2|) because sizes are different


#install.packages("BSDA")
library(BSDA)
z.test(my_transform(pcc_pwy),mu=mean(H0),sigma.x=sd(H0)) #p-value < 2.2e-16
z.test(my_transform(pcc_pcom),mu=mean(H0),sigma.x=sd(H0)) #p-value < 2.2e-16
z.test(my_transform(pcc_pwyANDpcom),mu=mean(H0),sigma.x=sd(H0)) #p-value < 2.2e-16
z.test(my_transform(pcc_5Annot),mu=mean(H0),sigma.x=sd(H0)) #p-value < 2.2e-16
