#Part I: Regulon
regulonAnnotation=list()

##pcc
start.time=Sys.time()


regulonAnnotation[1]=dist_TF_cumsum(data=All_Data_NAimputed,attribute_list = idInNichols_regulatorID,dist_metric=pcc_dist)
end.time=Sys.time()
end.time-start.time #Time difference of 10.3136 mins

##spearman
start.time=Sys.time()

regulonAnnotation[2]=dist_TF_cumsum(data=All_Data_NAimputed,attribute_list = idInNichols_regulatorID,dist_metric=spearman_dist)
end.time=Sys.time()
end.time-start.time 

##MHD3
start.time=Sys.time()
regulonAnnotation[3]=dist_TF_cumsum(data=Ternary_Data_324cutff_NAremoved,attribute_list = idInNichols_regulatorID,dist_metric=modified_hamming_distance_3)
end.time=Sys.time()
end.time-start.time #Time difference of 10.59211 mins (run on my PC)


##MHD3 on collapsed conds
start.time=Sys.time()
regulonAnnotation[4]=dist_TF_cumsum(data=Ternary_Data_324cutff_condCollapsed,attribute_list = idInNichols_regulatorID,dist_metric=modified_hamming_distance_3)
end.time=Sys.time()
end.time-start.time #Time difference of 11.05224 mins

##(Might wanna do this in the future)mutual info on ternary

##(Might wanna do this in the future)mutual info on collapsed conds



#save(regulonAnnotation,file="Data/sourced/regulonAnnotation.RData")



#Part II: Operon
operonAnnotation=list()

start.time=Sys.time()

pcc_dist=function(data){
  pairwise_cor=cor(t(data)) #Have to transpose to get correlations of the rows
  pcc_dist=1-abs(pairwise_cor) #I use abs() because anti-correlation is better than not having any relationships
  pcc_dist
}

operonAnnotation[1]=dist_TF_cumsum(data=All_Data_NAimputed,attribute_list = idInNichols_operonID,dist_metric=pcc_dist)
end.time=Sys.time()
end.time-start.time



start.time=Sys.time()

spearman_dist=function(data){
  if(!require(factoextra)){
    install.packages("factoextra")
    library(factoextra)
  }
  
  spearman=1-abs(cor(t(data),method="spearman",use="pairwise.complete.obs"))
  spearman
}

operonAnnotation[2]=dist_TF_cumsum(data=All_Data_NAimputed,attribute_list = idInNichols_operonID,dist_metric=spearman_dist)
end.time=Sys.time()
end.time-start.time


start.time=Sys.time()
operonAnnotation[3]=dist_TF_cumsum(data=Ternary_Data_324cutff_NAremoved,attribute_list = idInNichols_operonID,dist_metric=modified_hamming_distance_3)
end.time=Sys.time()
end.time-start.time 

start.time=Sys.time()
operonAnnotation[4]=dist_TF_cumsum(data=Ternary_Data_324cutff_condCollapsed,attribute_list = idInNichols_operonID,dist_metric=modified_hamming_distance_3)
end.time=Sys.time()
end.time-start.time

##(Might wanna do this in the future)mutual info on ternary

##(Might wanna do this in the future)mutual info on collapsed conds

#save(operonAnnotation,file="Data/sourced/operonAnnotation.RData")
