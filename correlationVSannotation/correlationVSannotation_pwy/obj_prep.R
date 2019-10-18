#Prepare objs that take some long time to process before graphing 

# Create a list: id - all pwy annotations

pwy.annotations=list()
for(i in 1:3979){
  TFvector=(ECK.Pathway$ids==i)
  pwy.annotations[[i]]=ECK.Pathway$`Data`[TFvector]
}

# Get pairwise ids by combn()
strain1strain2=as.data.frame(t(combn(1:3979,2)))
names(strain1strain2)=c("strain1","strain2")

# Create the final object
start=Sys.time()
bi.pwy.annot=c()
for(i in 1:dim(strain1strain2)[1]){
  strain1=strain1strain2[i,1]
  strain2=strain1strain2[i,2]
  bi.pwy.annot[i]=ifelse(sum(pwy.annotations[[strain1]] %in% pwy.annotations[[strain2]])>=1 & 
                           !(NA %in% pwy.annotations[[strain1]]) &
                           !(NA %in% pwy.annotations[[strain2]]) ,1,0)
}
end=Sys.time()
end-start #Time difference of 3.873857 mins

strain1strain2.samePWY_Annot=cbind(strain1strain2,bi.pwy.annot)
#save(strain1strain2.samePWY_Annot,file="Data/strain1strain2.samePWY_Annot.RData")



# Quantitative data preparation
start.time=Sys.time()
strain1strain2.samePWY_Annot.TF=merge(abs.sort.pcc.sql.NoIdent,strain1strain2.samePWY_Annot,by=c("strain1","strain2"))
end.time=Sys.time()
end.time-start.time #Time difference of 2.562144 mins
#Sort by pcc (decreasing)
strain1strain2.samePWY_Annot.TF=strain1strain2.samePWY_Annot.TF[order(strain1strain2.samePWY_Annot.TF$Pearson.Correlation.Coefficient,decreasing=T),]

pcc_TF=strain1strain2.samePWY_Annot.TF$bi.pwy.annot
pcc_coef=strain1strain2.samePWY_Annot.TF$Pearson.Correlation.Coefficient


# MHD (no punishment)
start.time=Sys.time()
strain1strain2.hammingNo0.samePWY_Annot=merge(binary_coef_hamming_No0,strain1strain2.samePWY_Annot,by=c("strain1","strain2"))
end.time=Sys.time()
end.time-start.time #Time difference of 1.002381 mins
#Sort by hamming distance (decreasing) .In hamming_No0 coefficient the larger the number is the more similar the 2 sequences are
#=> The real "distance" should be 324 (total elements compared) - binary_coef_hamming_No0
strain1strain2.hammingNo0.samePWY_Annot=strain1strain2.hammingNo0.samePWY_Annot[order(strain1strain2.hammingNo0.samePWY_Annot$'Hamming Distance',decreasing=T),]

MHDnoPunishment_TF=strain1strain2.hammingNo0.samePWY_Annot$bi.pwy.annot
MHDnoPunishment_coef=strain1strain2.hammingNo0.samePWY_Annot$'Hamming Distance'

# MHD (-1 for (1,-1) pairs)
start.time=Sys.time()
strain1strain2.hammingNo0_corrected.samePWY_Annot=merge(binary_coef_hamming_No0_corrected,strain1strain2.samePWY_Annot,by=c("strain1","strain2"))
end.time=Sys.time()
end.time-start.time #Time difference of 1.002381 mins
#Sort by hamming distance (decreasing) .In hamming_No0_corrected coefficient the larger the number is the more similar the 2 sequences are
#=> The real "distance" should be 324 (total elements compared) - binary_coef_hamming_No0_corrected
strain1strain2.hammingNo0_corrected.samePWY_Annot=strain1strain2.hammingNo0_corrected.samePWY_Annot[order(strain1strain2.hammingNo0_corrected.samePWY_Annot$'Hamming Distance_no0_corrected',decreasing=T),]
#save(strain1strain2.hammingNo0_corrected.samePWY_Annot,file="Data/sourced/strain1strain2.hammingNo0_corrected.samePWY_Annot.RData")

MHD_TF=strain1strain2.hammingNo0_corrected.samePWY_Annot$bi.pwy.annot
MHD_coef=strain1strain2.hammingNo0_corrected.samePWY_Annot$'Hamming Distance_no0_corrected'


#only genes in Nichols' figS1 pwy

##quantitative
inPwy.pcc_TF=pwy.pcc$`At least 1 same annotation`
inPwy.pcc_coef=pwy.pcc$Pearson.Correlation.Coefficient

##qualitative

## MHD (no punishment)

##get the table by merging 
inpwys_strain1strain2.hammingNo0.samePWY_Annot=merge(pwy.pcc[,c(1:2,4)],binary_coef_hamming_No0,by=c("strain1","strain2"))
#Sort by hamming distance (decreasing) .In hamming_No0 coefficient the larger the number is the more similar the 2 sequences are
inpwys_strain1strain2.hammingNo0.samePWY_Annot=inpwys_strain1strain2.hammingNo0.samePWY_Annot[order(inpwys_strain1strain2.hammingNo0.samePWY_Annot$'Hamming Distance',decreasing=T),]

##(!) Solve this: ids with the same genotypes will cause problem in this analysis (Eg. duplicated strains)

inPwy.MHDnoPunishment_TF=inpwys_strain1strain2.hammingNo0.samePWY_Annot$`At least 1 same annotation`
inPwy.MHDnoPunishment_coef=inpwys_strain1strain2.hammingNo0.samePWY_Annot$'Hamming Distance'


# MHD (-1 for (1,-1) pairs)
##get the table by merging 
inpwys_strain1strain2.hammingNo0_corrected.samePWY_Annot=merge(pwy.pcc[,c(1:2,4)],binary_coef_hamming_No0_corrected,by=c("strain1","strain2"))
#Sort by hamming distance (decreasing) .In hamming_No0 coefficient the larger the number is the more similar the 2 sequences are
inpwys_strain1strain2.hammingNo0_corrected.samePWY_Annot=inpwys_strain1strain2.hammingNo0_corrected.samePWY_Annot[order(inpwys_strain1strain2.hammingNo0_corrected.samePWY_Annot$'Hamming Distance',decreasing=T),]

##(!) Solve this: ids with the same genotypes will cause problem in this analysis (Eg. duplicated strains)

inPwy.MHD_TF=inpwys_strain1strain2.hammingNo0_corrected.samePWY_Annot$`At least 1 same annotation`
inPwy.MHD_coef=inpwys_strain1strain2.hammingNo0_corrected.samePWY_Annot$'Hamming Distance'



save(pcc_TF,pcc_coef,
     MHDnoPunishment_TF,
     MHDnoPunishment_coef,
     MHD_TF,
     MHD_coef,
     inPwy.pcc_TF,
     inPwy.pcc_coef,
     inPwy.MHDnoPunishment_TF,
     inPwy.MHDnoPunishment_coef,
     inPwy.MHD_TF,
     inPwy.MHD_coef,
     file="Data/sourced/pwy_TFsAndCoefs.RData")



# For distances using NA-imputed data

##All genes using Nichols' fig.S1 pwys
dat=All_Data_NAimputed
attribute_list=attr_list(ECK.Pathway$ids,ECK.Pathway$Data)

start.time = Sys.time()
pwy_NAimputed.listOfResults=dist_TF_cumsum_pcc_MHD(dat,ternary = T,thresh = 3.463,attribute_list)
end.time = Sys.time()
end.time - start.time #Time difference of 12.93831 mins


start.time = Sys.time()
pwy_NAimputed.listOfResults[3]=dist_TF_cumsum_spearman(dat,attribute_list)
end.time = Sys.time()
end.time - start.time #Time difference of 14.73373 mins


##This is for a test MHD:
dat=All_Data_NAimputed
attribute_list=attr_list(ECK.Pathway$ids,ECK.Pathway$Data)

start.time = Sys.time()
pwy_NAimputed.listOfResults[4]=dist_TF_cumsum_MHD4(dat,ternary = T,thresh = 3.463,attribute_list)
end.time = Sys.time()
end.time - start.time  #Time difference of 14.13567 mins



##This is for a test MHD:
dat=All_Data_NAimputed
attribute_list=attr_list(ECK.Pathway$ids,ECK.Pathway$Data)

start.time = Sys.time()
pwy_NAimputed.listOfResults[5]=dist_TF_cumsum_MHDs(dat,modified_hamming_distance_5,ternary = T,thresh = 3.463,attribute_list)
colnames(pwy_NAimputed.listOfResults[[5]])[3]="MHD5"
end.time = Sys.time()
end.time - start.time  #Time difference of 13.70093 mins
#save(pwy_NAimputed.listOfResults,file="Data/sourced/pwy_NAimputed.listOfResults.RData")



##Only genes in Nichols' fig.S1 pwys
strainInPwy=as.character(unique(c(pwy.pcc$strain1,pwy.pcc$strain2)))
dat=All_Data_NAimputed[strainInPwy,]
dim(dat)
index=(ECK.Pathway$ids %in% strainInPwy)
attribute_list=attr_list(ECK.Pathway$ids[index],ECK.Pathway$Data[index])

start.time = Sys.time()
inpwy_NAimputed.listOfResults=dist_TF_cumsum_pcc_MHD(dat,ternary = T,thresh = 3.463,attribute_list)
end.time = Sys.time()
end.time - start.time #Time difference of 24.24847 secs
#save(inpwy_NAimputed.listOfResults,file="Data/sourced/inpwy_NAimputed.listOfResults.RData")


























