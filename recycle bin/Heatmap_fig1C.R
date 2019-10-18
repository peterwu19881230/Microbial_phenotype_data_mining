#Goal: Do unit testing of Heatmap() function using part of Nichols' data


#Load the required libraries
library(ComplexHeatmap)
library(circlize)





#Import the data from .csv and tell the read.csv() that 1st column, 1st row are there for the names of the data
New_Data<-read.csv("allData_rows_fixed.csv",row.names = 1,header = TRUE)



#Transpose for cor() to work
#t_New_Data=t(New_Data)


#Subset
#(I don't know why there is an "unspecified" tag after SDS and some other drugs, but those seem to be just simply SDS or that particular treatment)
#row names in the R matrix were automatically truncated
conditions<-c("BILE.0.5....UNSPECIFIED","BILE.1.0....UNSPECIFIED","BILE.2.0....UNSPECIFIED","DEOXYCHOLATE.2.0....UNSPECIFIED","DEOXYCHOLATE.0.5....UNSPECIFIED","DEOXYCHOLATE.0.1....UNSPECIFIED",
"SDS0.5..EDTA0.5..","SDS1.0..EDTA0.5..","SDS0.5..EDTA0.1..",
"SDS.0.5....UNSPECIFIED","SDS.1.0....UNSPECIFIED","SDS.2.0....UNSPECIFIED","SDS.3.0....UNSPECIFIED","SDS.4.0....UNSPECIFIED",
"BENZALKONIUM.10...UNSPECIFIED","BENZALKONIUM.25...UNSPECIFIED",
"DIBUCAINE.0.4...UNSPECIFIED","DIBUCAINE.0.8...UNSPECIFIED","DIBUCAINE.1.2...UNSPECIFIED",
"NOVOBIOCIN.30...UNSPECIFIED",
"TRICLOSAN.0.05...UNSPECIFIED",
"ACRIFLAVINE.10...UNSPECIFIED",
"ETHIDIUMBROMIDE.50...UNSPECIFIED",
"ACRIFLAVINE.2...UNSPECIFIED",
"ETHIDIUMBROMIDE.10...UNSPECIFIED","ETHIDIUMBROMIDE.2...UNSPECIFIED",
"PROPIDIUMIODIDE.20...UNSPECIFIED","PROPIDIUMIODIDE.50...UNSPECIFIED",
"MINOCYCLINE.0.2...UNSPECIFIED","MINOCYCLINE.0.5...UNSPECIFIED","MINOCYCLINE.1.0...UNSPECIFIED",
"PUROMYCIN.25...UNSPECIFIED","PUROMYCIN.5...UNSPECIFIED",
"PYOCYANIN.10.0...UNSPECIFIED",
"MITOMYCINC.0.1...UNSPECIFIED",
"STREPTONIGRIN.0.5...UNSPECIFIED"
)
#strains<-c("ECK3622-RFAQ","ECK3617-RFAI","ECK3618-RFAB","ECK3620-RFAP","ECK3621-RFAG","ECK3610-RFAF","ECK3609-RFAD","ECK0223-LPCA","ECK3042-RFAE" )
strains<-c("ECK2272-NUOL","ECK2281-NUOB","ECK2273-NUOK","ECK2282-NUOA","ECK2274-NUOJ","ECK2276-NUOH","ECK2270-NUON","ECK2278-NUOF","ECK2271-NUOM","ECK2279-NUOE")
New_Data<-New_Data[strains,conditions]


#This trys to mimic the orginal settings in Nichols' as much as possible
Heatmap(New_Data,col=colorRamp2(c(-10,0,10),c("green","black","red")),cluster_rows=FALSE,cluster_columns=FALSE)
#Heatmap(t_New_Data,col=colorRamp2(c(-10,0,10),c("green","black","red"),cluster_rows=FALSE,cluster_columns=FALSE)) <-No idea why this didn't work. It worked in one case I tested before





