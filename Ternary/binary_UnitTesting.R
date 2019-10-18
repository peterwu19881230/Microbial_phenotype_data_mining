#Goal: Convert all data into binary (-1,0,1) using R

#Load the required libraries
#(Have to do this first: 1. source("https://bioconductor.org/biocLite.R") 2. biocLite("ComplexHeatmap"))
library(ComplexHeatmap)
library(circlize)

#Import the data from .csv and tell the read.csv() that 1st column, 1st row are there for the names of the data
All_Data<-read.csv("allData_rows_fixed.csv",row.names = 1,header = TRUE)


#subset to test first
test_original<-All_Data[1:500,1:324]
test<-All_Data[1:500,1:324]

test

#This is not the most elegant way...(and the conversion to 0s should be done first)
test[test>-3&test<3]<-0
test[test<=-3]<--1
test[test>=3]<-1


Heatmap(test)
Heatmap(test_original)
