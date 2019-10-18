#Pairwise correlations for distances used for different types of data 
#(!)(if different types of  quantitative/qualitative data are used, they are not that different: slightly different cutoffs or whether NAs are removed)



listOfResults=pwy_NAimputed.listOfResults
pcc=pwy_NAimputed.listOfResults[[1]]$pcc
mhd=pwy_NAimputed.listOfResults[[2]]$MHD3
mhd_collapsedConditions=pwy_collapsedCond.MHD3[[1]]$MHD3
load("Data/mutual_info_ternaryData_table.RData")
mutual_info=mutual_info_ternaryData_table[,3]

cor(pcc,mhd,method = "spearman") #0.2635486
cor(pcc,mhd_collapsedConditions,method = "spearman") #0.3810959
cor(pcc,mutual_info,method = "spearman") #0.010931
cor(mhd,mhd_collapsedConditions,method = "spearman") #0.6857363
cor(mhd,mutual_info,method = "spearman") #-0.041552
cor(mhd_collapsedConditions,mutual_info,method = "spearman") #-0.02875783
