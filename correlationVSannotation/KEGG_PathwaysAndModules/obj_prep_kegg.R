#Goal: obj preparation for KEGG module analysis

# 1. Get the obj for quantitative data
library(dplyr)
load("Data/KEGGmodulesForNichols.RData")

start.time = Sys.time()
attribute_list=KEGGmodulesForNichols
listOfResults=dist_TF_cumsum_old(All_Data_NAimputed[names(attribute_list),],attribute_list)
end.time = Sys.time()
end.time - start.time #Time difference of 25.4773 mins
kegg_quantitative=listOfResults
#save(kegg_quantitative,file="Data/sourced/kegg_quantitative.RData")

## I prepare this TF vector to subset the strains annotated with modules
inModulesOrNot=sapply(KEGGmodulesForNichols,FUN=function(strain){
  ifelse(identical(strain,character(0)),0,1)  
})

## For strains annotated with modules
attribute_list=KEGGmodulesForNichols[as.logical(inModulesOrNot)] ## If as.logical is not used, it will subset the 1st element for multiple times. 0 won't get anything
listOfResults=dist_TF_cumsum_old(All_Data_NAimputed[names(attribute_list),],attribute_list)
keggInModules_quantitative=listOfResults
#save(keggInModules_quantitative,file="Data/sourced/keggInModules_quantitative.RData")

# 2. Get the objs for ternary data

##get MHD distance obj
#MHD_ternary ##(-1,-1) (1,1) +1
#MHD2_ternary ##(-1,-1) (1,1) +1 #(-1,1) (1,-1) -1 #(0,-1) (-1,0) (0,1) (1,0) -0.5
#MHD3_ternary ##(-1,-1) (1,1) +1 #(-1,1) (1,-1) -1

#melt the distance
MHDs_prosessing=function(anMHD){
  MHD_melted=melt_dist(as.dist(anMHD))
  
  #Bind the TFvector (this part can possibly re-written by parSapply)
  distance_table=MHD_melted
  sameORnot=c()
  for(i in 1:dim(distance_table)[1]){
    id1=distance_table[i,1]
    id2=distance_table[i,2]
    
    sameORnot[i]=ifelse( sum( attribute_list[[id1]] %in% attribute_list[[id2]])>=1 & #The way attribute_list[[id1]] works is like attribute_list$id1, where id1 is a string. I haven't tested the function when id1 is a number 
                           !(anyNA(attribute_list[[id1]])) & 
                           !(anyNA(attribute_list[[id2]])),
                         1,0)
  }
  
  distance_table_TF=cbind(distance_table,sameORnot)
  distance_table_TF=distance_table_TF[order(distance_table_TF$Distance),]
  
  return(distance_table_TF)
}

start.time = Sys.time()
MHDdistance_table=MHDs_prosessing(MHD_ternary)
MHD2distance_table=MHDs_prosessing(MHD2_ternary)
MHD3distance_table=MHDs_prosessing(MHD3_ternary)
cumSum_MHDdistance_table_TF=cumsum(MHDdistance_table$sameORnot)
cumSum_MHD2distance_table_TF=cumsum(MHD2distance_table$sameORnot)
cumSum_MHD3distance_table_TF=cumsum(MHD3distance_table$sameORnot)
end.time = Sys.time()
end.time - start.time #Time difference of 26.73598 mins

## For strains annotated with modules
strains=names(inModulesOrNot)[as.logical(inModulesOrNot)]
MHDdistanceInModule_table=MHDs_prosessing(MHD_ternary[strains,strains])
MHD2distanceInModule_table=MHDs_prosessing(MHD2_ternary[strains,strains])
MHD3distanceInModule_table=MHDs_prosessing(MHD3_ternary[strains,strains])
cumSum_MHDdistanceInModule_table_TF=cumsum(MHDdistanceInModule$sameORnot)
cumSum_MHD2distanceInModule_table_TF=cumsum(MHD2distanceInModule$sameORnot)
cumSum_MHD3distanceInModule_table_TF=cumsum(MHD3distanceInModule$sameORnot)


#(Note: if saving too may obj in 1 file, there will not be enough memory to load the file)
save(kegg_quantitative,file="Data/sourced/kegg_quantitative.RData")
save(keggInModules_quantitative,file="Data/sourced/keggInModules_quantitative.RData")
save(MHDdistance_table,file="Data/sourced/MHDdistance_table.RData")
save(MHD2distance_table,file="Data/sourced/MHD2distance_table.RData")
save(MHD3distance_table,file="Data/sourced/MHD3distance_table.RData")

#save(cumSum_MHDdistance_table_TF,cumSum_MHD2distance_table_TF, cumSum_MHD3distance_table_TF,cumSum_MHDdistanceInModule_table_TF,cumSum_MHD2distanceInModule_table_TF,cumSum_MHD3distanceInModule_table_TF,
#     file="Data/sourced/obj_prep_kegg_2.RData")



#When data = using only 114 treatments collapsed by 324 conditions:
start.time = Sys.time()
attribute_list=KEGGmodulesForNichols
listOfResults=dist_TF_cumsum_MHD3(Ternary_Data_condCollapsed[names(attribute_list),],attribute_list=attribute_list)
end.time = Sys.time()
end.time - start.time #Time difference of 21.22536 mins
kegg_Ternary_Data_condCollapsed_MHD3=listOfResults
#save(kegg_Ternary_Data_condCollapsed_MHD3,file="Data/sourced/kegg_Ternary_Data_condCollapsed_MHD3.RData")







