#Goal:
#1.Heatmap: reproduce 4A (using similar colors)
#2.Heatmap: use correlation coefficients calculated from binary.csv -> doesn't work at all! All the values are NAs (Heatmap requres at least 2 distinct values)




#Create a list to subset
fig4A.list=GetIndices(c("ECK0457-ACRA","ECK0456-ACRB","ECK1524-MARA","ECK1525-MARB","ECK1523-MARR"),ECK_1st_table$sorted_ECK_missing_gene_names_added)
##GetIndices is a self-defined function

fig4A.matrix=cor_strains[fig4A.list,fig4A.list]

#Load the required libraries
library(ComplexHeatmap)
library(circlize)


pdf('fig4A.pdf')
#Result of 1
Heatmap(fig4A.matrix,col=colorRamp2(c(-1,0,1),c("cyan","black","yellow")),cluster_rows=FALSE, cluster_columns=FALSE )
dev.off()


