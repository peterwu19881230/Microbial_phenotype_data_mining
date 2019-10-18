#Goal: Do Heatmap() function using Nichols' data->binary and output it to a pdf file


#Load the required libraries
library(ComplexHeatmap)
library(circlize)





#Import the data from .csv and tell the read.csv() that 1st column, 1st row are there for the names of the data
New_Data<-read.csv("Data/binary_Nichols.csv",row.names = 1,header = TRUE)


#col_names and row_names are not easily seen on the plot (too dense), so I remove them
rownames(New_Data)<-c()
colnames(New_Data)<-c()


#This trys to mimic the orginal settings in Nichols' as much as possible
#However, I don't know what's the clustering method in detail
pdf('nocluster_binary_Nichols_heatmap.pdf')

#I don't know why this causes memory problem Even the original one (with original phenotype scores) doesn't
#Heatmap(t_New_Data,col=colorRamp2(c(-1,0,1),c("green","black","red")))
Heatmap(New_Data,col=colorRamp2(c(-1,0,1),c("blue","white","red")),cluster_rows=FALSE,cluster_columns=FALSE)


#A small portion of it:
pdf('small_nocluster_binary_Nichols_heatmap.pdf')
Heatmap(New_Data[1:50,1:50],col=colorRamp2(c(-1,0,1),c("blue","white","red")),cluster_rows=FALSE,cluster_columns=FALSE)


