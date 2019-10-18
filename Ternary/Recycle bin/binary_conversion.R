#Goal: Convert all data into binary (-1,0,1) using R


#This is not the most elegant way
matrix<-matrix(,3979,324)
rownames(matrix)<-rownames(All_Data)
colnames(matrix)<-colnames(All_Data)
matrix[All_Data>-3.463&All_Data<3.463]<-0
matrix[All_Data<=-3.463]<--1
matrix[All_Data>=3.463]<-1

###
#Another way to convert original data to binary...(the conversion to 0s should be done first):
#matrix<-All_Data
#thres<-3.463 #Set the threshold to 3.463
#matrix[matrix>(thres*(-1))&matrix<thres]<-0
#matrix[matrix<(thres*(-1))]<--1
#matrix[matrix>thres]<-1

#Test if the converted matrix is the same as the original Binary_Data:
#Binary_Data[1:10,1:10]==matrix[1:10,1:10]
###



write.csv(matrix, file = "Data/binary_Nichols.csv")

#Load the required libraries for heatmaps
#library(ComplexHeatmap)
#library(circlize)

# (On my computer) Running Heatmap() with the part that converts values to binary will return an error: Error: node stack overflow
#Heatmap(All_Data)
#Heatmap(All_Data_original)
