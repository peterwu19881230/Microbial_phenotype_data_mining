#Goal:
#1.Heatmap: reproduce 5B (using similar colors)
#2.Heatmap: replace scores with >3: 1, -3<socre<3: 0, <3: -1





#Find what's in the fig5B
load("Data/ECK_1st_table.RData")
columns<-GetIndices(c("ECK2898-GCVP","ECK2900-GCVT","ECK2899-GCVH","ECK2548-GLYA"),ECK_1st_table$sorted_ECK_missing_gene_names_added)  

#I don't know why R has to covert the row name to this weired format (all unacceptable characters are all converted to ".")
rows<-c("SULFAMETHIZOLE.100...UNSPECIFIED","SULFAMETHIZOLE.200...UNSPECIFIED","SULFAMETHIZOLE.300...UNSPECIFIED","SULFAMONOMETHOXINE.50...UNSPECIFIED","SULFAMONOMETHOXINE.100...UNSPECIFIED","TRIMETHOPRIM.0.1.SULFAMETHIZOLE.50...UNSPECIFIED","TRIMETHOPRIM.0.1...UNSPECIFIED","TRIMETHOPRIM.0.2...UNSPECIFIED","TRIMETHOPRIM.0.3...UNSPECIFIED","TRIMETHOPRIM.0.4...UNSPECIFIED")


#Subset
fig5B.matrix<-All_Data[columns,rows]


#Conversion for result2
binary<-fig5B.matrix
binary[binary>-3&binary<3]<-0
binary[binary<=-3]<--1
binary[binary>=3]<-1


#Info:

#sulfa: sulfanimides (In the paper it shows 5 conditions of 1 single compound. However I found 2 in the excel)
#SULFAMETHIZOLE-100 - UNSPECIFIED
#SULFAMETHIZOLE-200 - UNSPECIFIED
#SULFAMETHIZOLE-300 - UNSPECIFIED
#SULFAMONOMETHOXINE-50 - UNSPECIFIED
#SULFAMONOMETHOXINE-100 - UNSPECIFIED

#TMP+sulfa (In the paper 3 conditions were shown. However in the excel file only this condition was found)
#TRIMETHOPRIM-0.1,SULFAMETHIZOLE-50 - UNSPECIFIED

#TMP
#TRIMETHOPRIM-0.1 - UNSPECIFIED
#TRIMETHOPRIM-0.2 - UNSPECIFIED
#TRIMETHOPRIM-0.3 - UNSPECIFIED
#TRIMETHOPRIM-0.4 - UNSPECIFIED




#Load the required libraries
library(ComplexHeatmap)
library(circlize)


pdf('fig5B.pdf') #pdf() and jpeg() won't give the proper labeling (text descriptions go out of the pages)
#Result of 1
Heatmap(fig5B.matrix,col=colorRamp2(c(-8,0,8),c("green","black","red")),cluster_rows=FALSE, cluster_columns=FALSE, column_names_gp = gpar(fontsize = 6) )

#Result of 2
Heatmap(binary,col=colorRamp2(c(-1,0,1),c("green","black","red")),cluster_rows=FALSE, cluster_columns=FALSE, column_names_gp = gpar(fontsize = 6)  )
dev.off()


#cluster based only on the drugs in fig5B
library(dplyr)
All_Data_fig5B<-All_Data[,c("SULFAMETHIZOLE.100...UNSPECIFIED","SULFAMETHIZOLE.200...UNSPECIFIED","SULFAMETHIZOLE.300...UNSPECIFIED","SULFAMONOMETHOXINE.50...UNSPECIFIED","SULFAMONOMETHOXINE.100...UNSPECIFIED","TRIMETHOPRIM.0.1.SULFAMETHIZOLE.50...UNSPECIFIED","TRIMETHOPRIM.0.1...UNSPECIFIED","TRIMETHOPRIM.0.2...UNSPECIFIED","TRIMETHOPRIM.0.3...UNSPECIFIED","TRIMETHOPRIM.0.4...UNSPECIFIED")]
All_Data_fig5B %>% pcc_dist %>% hclust %>% plot



