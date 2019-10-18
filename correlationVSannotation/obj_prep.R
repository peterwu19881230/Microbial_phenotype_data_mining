#Goal: for obj preparation
# Quantitative data preparation
##pwy
start.time=Sys.time()
strain1strain2.pcc.pwyTF=merge(abs.sort.pcc.sql.NoIdent,strain1strain2.samePWY_Annot,by=c("strain1","strain2"))
end.time=Sys.time()
end.time-start.time #Time difference of 2.562144 mins

##pcomplex
##take abs of pcc
strain1strain2.pcc.pcomTF=strain1strain2.pcc.samePTcomplex_Annot
strain1strain2.pcc.pcomTF$Pearson.Correlation.Coefficient=abs(strain1strain2.pcc.pcomTF$Pearson.Correlation.Coefficient)


##pwy+pcomplex
start.time=Sys.time()
strain1strain2.pcc.samePwyAnnot.samePcomAnnot=merge(strain1strain2.pcc.pwyTF,strain1strain2.pcc.pcomTF,by=c("strain1","strain2"))
end.time=Sys.time()
end.time-start.time #Time difference of 1.128728 mins


pcc_TF=apply(strain1strain2.pcc.samePwyAnnot.samePcomAnnot[,c("bi.pwy.annot","bi.ptcomplex.annot")],1,FUN=function(row){
  if(row[1]==0 & row[2]==0){
    return(0)
  }else{1}
})

##A table that is: strain1 - strain2 - pcc - PwyPcom_TF
if(!exists("pccPwyPcom")){
pccPwyPcom=cbind(strain1strain2.pcc.samePwyAnnot.samePcomAnnot[,1:3],pcc_TF)
pccPwyPcom=pccPwyPcom[order(pccPwyPcom$Pearson.Correlation.Coefficient.x,decreasing=T),]##sort by pcc
save(pccPwyPcom,file="Data/pccPwyPcom.RData")
}

if(!exists("hammingPwyPcom")){
  # Binary data preparation
  ##pwy
  start.time=Sys.time()
  strain1strain2.hammingNo0.samePWY_Annot=merge(binary_coef_hamming_No0,strain1strain2.samePWY_Annot,by=c("strain1","strain2"))
  end.time=Sys.time()
  end.time-start.time #Time difference of 1.25747 mins
  
  ##pcomplex
  start.time=Sys.time()
  strain1strain2.hammingNo0.samePTcomplex_Annot=merge(binary_coef_hamming_No0,strain1strain2.pcc.samePTcomplex_Annot[,-3],by=c("strain1","strain2"))
  end.time=Sys.time()
  end.time-start.time #Time difference of 1.246467 mins
  
  ##pwy+pcomplex
  start.time=Sys.time()
  strain1strain2.hammingNo0.samePwyAnnot.samePcomAnnot=merge(strain1strain2.hammingNo0.samePWY_Annot,strain1strain2.hammingNo0.samePTcomplex_Annot,by=c("strain1","strain2"))
  end.time=Sys.time()
  end.time-start.time #Time difference of 42.60484 secs
  
  hammingNo0_TF=apply(strain1strain2.hammingNo0.samePwyAnnot.samePcomAnnot[,c("bi.pwy.annot","bi.ptcomplex.annot")],1,FUN=function(row){
    if(row[1]==0 & row[2]==0){
      return(0)
    }else{1}
  })
  
  
  ##A table that is: strain1 - strain2 - hamNo0 - PwyPcom_TF
  hammingPwyPcom=cbind(strain1strain2.hammingNo0.samePwyAnnot.samePcomAnnot[,1:3],hammingNo0_TF)
  hammingPwyPcom=hammingPwyPcom[order(hammingPwyPcom$`Hamming Distance.x`,decreasing=T),]##sort by hamNo0
  save(hammingPwyPcom,file="Data/hammingPwyPcom.RData")
}


#Other distance methods to compare with PCC based + quantitative
strain1strain2.samePwyPtcomplex_Annot=pccPwyPcom[,c(1,2,4)] #strain1-strain2-(Same pwy/ptcom annotations or not)

'
start.time=Sys.time()

strain1strain2.distances.PwyPtcomTF=Reduce(function(x, y) merge(x, y, by=c("strain1","strain2")), 
list(strain1strain2.samePwyPtcomplex_Annot, 
sort.dist.All_Data_NAimputed.euclidean, 
sort.dist.All_Data_NAimputed.maximum,
sort.dist.All_Data_NAimputed.manhattan,
sort.dist.All_Data_NAimputed.canberra,
sort.dist.All_Data_NAimputed.binary,
sort.dist.All_Data_NAimputed.minkowski,
sort.dist.All_Data_NAimputed.spearman
))

end.time=Sys.time()
end.time-start.time #Time difference of 8.670135 mins

##reorder the columns and correct the name of the TF vector
strain1strain2.distances.PwyPtcomTF=strain1strain2.distances.PwyPtcomTF[,c(1,2,4,5,6,7,8,9,10,3)]
names(strain1strain2.distances.PwyPtcomTF)[10]="TF"

save(strain1strain2.distances.PwyPtcomTF,file="Data/strain1strain2.distances.PwyPtcomTF.RData")
'

#Different MHD (computed after NAs were imputed with 0)
'
strain1strain2.MHD=meltANDsort_dist(MHD_ternary)
names(strain1strain2.MHD)=c("strain1","strain2","MHD")
strain1strain2.MHD.samePwyPtcomplex_Annot=merge(strain1strain2.MHD,strain1strain2.samePwyPtcomplex_Annot,by=c("strain1","strain2"))
strain1strain2.MHD2=meltANDsort_dist(MHD2_ternary)
names(strain1strain2.MHD2)=c("strain1","strain2","MHD2")
strain1strain2.MHD2.samePwyPtcomplex_Annot=merge(strain1strain2.MHD2,strain1strain2.samePwyPtcomplex_Annot,by=c("strain1","strain2"))
strain1strain2.MHD3=meltANDsort_dist(MHD3_ternary)
names(strain1strain2.MHD3)=c("strain1","strain2","MHD3")
strain1strain2.MHD3.samePwyPtcomplex_Annot=merge(strain1strain2.MHD3,strain1strain2.samePwyPtcomplex_Annot,by=c("strain1","strain2"))

##correct a colname
names(sort.strain1strain2.MHD.samePwyPtcomplex_Annot)[4]="TF" 
names(sort.strain1strain2.MHD2.samePwyPtcomplex_Annot)[4]="TF" 
names(sort.strain1strain2.MHD3.samePwyPtcomplex_Annot)[4]="TF" 


save(strain1strain2.MHD.samePwyPtcomplex_Annot,
strain1strain2.MHD2.samePwyPtcomplex_Annot,
strain1strain2.MHD3.samePwyPtcomplex_Annot,
file="Data/strain1strain2.MHDs.samePwyPtcomplex_Annot.RData")
'


##For only genes in pathways
##pwy.pcc is used
'
start.time=Sys.time()

strain1strain2.distances.inPwyTF=Reduce(function(x, y) merge(x, y, by=c("strain1","strain2")), 
list(pwy.pcc[,c(1,2,4)], 
sort.dist.All_Data_NAimputed.euclidean, 
sort.dist.All_Data_NAimputed.maximum,
sort.dist.All_Data_NAimputed.manhattan,
sort.dist.All_Data_NAimputed.canberra,
sort.dist.All_Data_NAimputed.binary,
sort.dist.All_Data_NAimputed.minkowski,
sort.dist.All_Data_NAimputed.spearman
))

end.time=Sys.time()
end.time-start.time Time difference of 1.63894 mins

##reorder the columns and correct the name of the TF vector
strain1strain2.distances.inPwyTF=strain1strain2.distances.inPwyTF[,c(1,2,4,5,6,7,8,9,10,3)]
names(strain1strain2.distances.inPwyTF)[10]="TF"

save(strain1strain2.distances.inPwyTF,file="Data/strain1strain2.distances.inPwyTF.RData")
'

##Different in-pathway MHD (computed after NAs were imputed with 0)
'
strain1strain2.MHD=meltANDsort_dist(MHD_ternary)
names(strain1strain2.MHD)=c("strain1","strain2","MHD")
pwy_strain1strain2.MHD.Annot=merge(strain1strain2.MHD,pwy.pcc[,c(1,2,4)],by=c("strain1","strain2"))
strain1strain2.MHD2=meltANDsort_dist(MHD2_ternary)
names(strain1strain2.MHD2)=c("strain1","strain2","MHD2")
pwy_strain1strain2.MHD2.Annot=merge(strain1strain2.MHD2,pwy.pcc[,c(1,2,4)],by=c("strain1","strain2"))
strain1strain2.MHD3=meltANDsort_dist(MHD3_ternary)
names(strain1strain2.MHD3)=c("strain1","strain2","MHD3")
pwy_strain1strain2.MHD3.Annot=merge(strain1strain2.MHD3,pwy.pcc[,c(1,2,4)],by=c("strain1","strain2"))


save(pwy_strain1strain2.MHD.Annot,
pwy_strain1strain2.MHD2.Annot,
pwy_strain1strain2.MHD3.Annot,
file="Data/pwy_strain1strain2.MHDs.Annot.RData")
'

#For only genes in protein complexes
'
start.time=Sys.time()

strain1strain2.distances.inPcomTF=Reduce(function(x, y) merge(x, y, by=c("strain1","strain2")), 
                                         list(inptcomplex_strain1strain2.pcc.samePTcomplex_Annot[,c(1,2,4)], 
                                              sort.dist.All_Data_NAimputed.euclidean, 
                                              sort.dist.All_Data_NAimputed.maximum,
                                              sort.dist.All_Data_NAimputed.manhattan,
                                              sort.dist.All_Data_NAimputed.canberra,
                                              sort.dist.All_Data_NAimputed.binary,
                                              sort.dist.All_Data_NAimputed.minkowski,
                                              sort.dist.All_Data_NAimputed.spearman
                                         ))

end.time=Sys.time()
end.time-start.time 

##reorder the columns and correct the name of the TF vector
strain1strain2.distances.inPcomTF=strain1strain2.distances.inPcomTF[,c(1,2,4,5,6,7,8,9,10,3)]
names(strain1strain2.distances.inPcomTF)[10]="TF"

save(strain1strain2.distances.inPcomTF,file="Data/strain1strain2.distances.inPcomTF.RData")
'

## Different in-ptcomplex MHD (computed after NAs were imputed with 0)
'
strain1strain2.MHD=meltANDsort_dist(MHD_ternary)
names(strain1strain2.MHD)=c("strain1","strain2","MHD")
pcom_strain1strain2.MHD.Annot=merge(strain1strain2.MHD,inptcomplex_strain1strain2.pcc.samePTcomplex_Annot[,c(1,2,4)],by=c("strain1","strain2"))

strain1strain2.MHD2=meltANDsort_dist(MHD2_ternary)
names(strain1strain2.MHD2)=c("strain1","strain2","MHD2")
pcom_strain1strain2.MHD2.Annot=merge(strain1strain2.MHD2,inptcomplex_strain1strain2.pcc.samePTcomplex_Annot[,c(1,2,4)],by=c("strain1","strain2"))

strain1strain2.MHD3=meltANDsort_dist(MHD3_ternary)
names(strain1strain2.MHD3)=c("strain1","strain2","MHD3")
pcom_strain1strain2.MHD3.Annot=merge(strain1strain2.MHD3,inptcomplex_strain1strain2.pcc.samePTcomplex_Annot[,c(1,2,4)],by=c("strain1","strain2"))


save(pcom_strain1strain2.MHD.Annot,
pcom_strain1strain2.MHD2.Annot,
pcom_strain1strain2.MHD3.Annot,
file="Data/pcom_strain1strain2.MHDs.Annot.RData")
'



#Load the objs being used
load("Data/strain1strain2.distances.PwyPtcomTF.RData")
load("Data/strain1strain2.MHDs.samePwyPtcomplex_Annot.RData")
load("Data/strain1strain2.distances.inPwyTF.RData")
load("Data/strain1strain2inPwy.MHD.samePwyAnnot.RData")
load("Data/pwy_strain1strain2.MHDs.Annot.RData")
load("Data/strain1strain2.distances.inPcomTF.RData")
load("Data/pcom_strain1strain2.MHDs.Annot.RData")

#For all gene pairs

TF=list() #sorted TF using different dists will be used (yet sum(TF) stays the same for every dist)
TF$pcc=pccPwyPcom$pcc_TF
TF$MHD=hammingPwyPcom$hammingNo0_TF[order(hammingPwyPcom$`Hamming Distance.x`)]
TF$euclidean=strain1strain2.distances.PwyPtcomTF$TF[order(strain1strain2.distances.PwyPtcomTF$euclidean)]
TF$maximum=strain1strain2.distances.PwyPtcomTF$TF[order(strain1strain2.distances.PwyPtcomTF$maximum)]
TF$manhattan=strain1strain2.distances.PwyPtcomTF$TF[order(strain1strain2.distances.PwyPtcomTF$manhattan)]
TF$canberra=strain1strain2.distances.PwyPtcomTF$TF[order(strain1strain2.distances.PwyPtcomTF$canberra)]
TF$binary=strain1strain2.distances.PwyPtcomTF$TF[order(strain1strain2.distances.PwyPtcomTF$binary)]
TF$minkowski=strain1strain2.distances.PwyPtcomTF$TF[order(strain1strain2.distances.PwyPtcomTF$minkowski)]
TF$spearman=strain1strain2.distances.PwyPtcomTF$TF[order(strain1strain2.distances.PwyPtcomTF$spearman)]

##Add different versions of MHD:
TF$MHD_noNA=strain1strain2.MHD.samePwyPtcomplex_Annot$TF[order(strain1strain2.MHD.samePwyPtcomplex_Annot$MHD)]
TF$MHD2_noNA=strain1strain2.MHD2.samePwyPtcomplex_Annot$TF[order(strain1strain2.MHD2.samePwyPtcomplex_Annot$MHD2)]
TF$MHD3_noNA=strain1strain2.MHD3.samePwyPtcomplex_Annot$TF[order(strain1strain2.MHD3.samePwyPtcomplex_Annot$MHD3)]


##Calculate cumulative sum for each distance
cumSum=lapply(TF,cumsum)

slopeForControl=sum(TF$pcc)/length(TF$pcc) #this slope stays the same for every other distances


#For gene pairs in pathways

##MHD computed with NAs
#strain1strain2inPwy.MHD.samePwyAnnot=merge(pwy.pcc[,c(1,2,4)],hammingPwyPcom[,1:3],by=c("strain1","strain2")) 
#save(strain1strain2inPwy.MHD.samePwyAnnot,file="Data/strain1strain2inPwy.MHD.samePwyAnnot.RData")


TF_pwy=list() #sorted TF using different dists will be used (yet sum(TF) stays the same for every dist)
TF_pwy$pcc=pwy.pcc$`At least 1 same annotation`
TF_pwy$MHD=strain1strain2inPwy.MHD.samePwyAnnot$`At least 1 same annotation`[order(strain1strain2inPwy.MHD.samePwyAnnot$`Hamming Distance.x`)]
TF_pwy$euclidean=strain1strain2.distances.inPwyTF$TF[order(strain1strain2.distances.inPwyTF$euclidean)]
TF_pwy$maximum=strain1strain2.distances.inPwyTF$TF[order(strain1strain2.distances.inPwyTF$maximum)]
TF_pwy$manhattan=strain1strain2.distances.inPwyTF$TF[order(strain1strain2.distances.inPwyTF$manhattan)]
TF_pwy$canberra=strain1strain2.distances.inPwyTF$TF[order(strain1strain2.distances.inPwyTF$canberra)]
TF_pwy$binary=strain1strain2.distances.inPwyTF$TF[order(strain1strain2.distances.inPwyTF$binary)]
TF_pwy$minkowski=strain1strain2.distances.inPwyTF$TF[order(strain1strain2.distances.inPwyTF$minkowski)]
TF_pwy$spearman=strain1strain2.distances.inPwyTF$TF[order(strain1strain2.distances.inPwyTF$spearman)]


##Add different versions of MHD:
TF_pwy$MHD_noNA=pwy_strain1strain2.MHD.Annot$`At least 1 same annotation`[order(pwy_strain1strain2.MHD.Annot$MHD)]
TF_pwy$MHD2_noNA=pwy_strain1strain2.MHD2.Annot$`At least 1 same annotation`[order(pwy_strain1strain2.MHD2.Annot$MHD2)]
TF_pwy$MHD3_noNA=pwy_strain1strain2.MHD3.Annot$`At least 1 same annotation`[order(pwy_strain1strain2.MHD3.Annot$MHD3)]


##Calculate cumulative sum for each distance
cumSum_pwy=lapply(TF_pwy,cumsum)

slopeForControl_inPwy=sum(TF_pwy$pcc)/length(TF_pwy$pcc) #this slope stays the same for every other distances


#For gene pairs in protein complexes


#MHD computed with NAs
#strain1strain2inPcom.MHD.samePcomAnnot=merge(inptcomplex_strain1strain2.pcc.samePTcomplex_Annot[,c(1,2,4)],hammingPwyPcom[,1:3],by=c("strain1","strain2")) 
#save(strain1strain2inPcom.MHD.samePcomAnnot,file="Data/strain1strain2inPcom.MHD.samePcomAnnot.RData")


TF_pcom=list() #sorted TF using different dists will be used (yet sum(TF) stays the same for every dist)
TF_pcom$pcc=inptcomplex_strain1strain2.pcc.samePTcomplex_Annot$bi.ptcomplex.annot[order(abs(inptcomplex_strain1strain2.pcc.samePTcomplex_Annot$Pearson.Correlation.Coefficient),decreasing = T)]
TF_pcom$MHD=strain1strain2inPcom.MHD.samePcomAnnot$bi.ptcomplex.annot[order(strain1strain2inPcom.MHD.samePcomAnnot$`Hamming Distance.x`)]
TF_pcom$euclidean=strain1strain2.distances.inPcomTF$TF[order(strain1strain2.distances.inPcomTF$euclidean)]
TF_pcom$maximum=strain1strain2.distances.inPcomTF$TF[order(strain1strain2.distances.inPcomTF$maximum)]
TF_pcom$manhattan=strain1strain2.distances.inPcomTF$TF[order(strain1strain2.distances.inPcomTF$manhattan)]
TF_pcom$canberra=strain1strain2.distances.inPcomTF$TF[order(strain1strain2.distances.inPcomTF$canberra)]
TF_pcom$binary=strain1strain2.distances.inPcomTF$TF[order(strain1strain2.distances.inPcomTF$binary)]
TF_pcom$minkowski=strain1strain2.distances.inPcomTF$TF[order(strain1strain2.distances.inPcomTF$minkowski)]
TF_pcom$spearman=strain1strain2.distances.inPcomTF$TF[order(strain1strain2.distances.inPcomTF$spearman)]


##Add different versions of MHD:
TF_pcom$MHD_noNA=pcom_strain1strain2.MHD.Annot$bi.ptcomplex.annot[order(pcom_strain1strain2.MHD.Annot$MHD)]
TF_pcom$MHD2_noNA=pcom_strain1strain2.MHD2.Annot$bi.ptcomplex.annot[order(pcom_strain1strain2.MHD2.Annot$MHD2)]
TF_pcom$MHD3_noNA=pcom_strain1strain2.MHD3.Annot$bi.ptcomplex.annot[order(pcom_strain1strain2.MHD3.Annot$MHD3)]


##Calculate cumulative sum for each distance
cumSum_pcom=lapply(TF_pcom,cumsum)

slopeForControl_inPcom=sum(TF_pcom$pcc)/length(TF_pcom$pcc) #this slope stays the same for every other distances



