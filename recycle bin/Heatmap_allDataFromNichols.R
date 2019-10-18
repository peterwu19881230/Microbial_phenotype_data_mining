#Goal: Do Heatmap() function using Nichols' data and output it to a pdf file

#Load the required libraries
library(ComplexHeatmap)
library(circlize)





#Import the data from .csv and tell the read.csv() that 1st column, 1st row are there for the names of the data
New_Data<-read.csv("allData_rows_fixed.csv",row.names = 1,header = TRUE)

#col_names and row_names are not easily seen on the plot (too dense), so I remove them
rownames(New_Data)<-c()
colnames(New_Data)<-c()




#This trys to mimic the orginal settings in Nichols' as much as possible
#However, I don't know what's the clustering method in detail
pdf('Nichols_heatmap.pdf')
Heatmap(New_Data)

jpeg('Nichols_heatmap.jpg')
Heatmap(New_Data)


#Show part of it
#!Actually this is wrong because it only shows a small portion of it
#pdf('small_Nichols_heatmap.pdf')
#Heatmap(New_Data[1:50,1:50],col=colorRamp2(c(-3,0,3),c("blue","white","red")))



