#Goal: Do unit testing of correlation coefficient using part of Nichols' data



#Import the data from .csv and tell the read.csv() that 1st column, 1st row are there for the names of the data
New_Data<-read.csv("SubData.csv",row.names = 1,header = TRUE)



#Transpose for cor() to work
t_New_Data=t(New_Data)

#Change i and j later to calculate c.c. of different strains. I deliberately picked 2 clomns with NA values
i<-"ECK2623-YFJK"
j<-"ECK1602-YDGC"
Sub_t_New_Data<-t_New_Data[,c(i,j)]


#Clean up NA values
good<-complete.cases(Sub_t_New_Data[,c(i,j)])
clean_Sub_t_New_Data<-Sub_t_New_Data[good,]


#Calculate the correlation coefficient
cor(clean_Sub_t_New_Data[,i],clean_Sub_t_New_Data[,j])


#plot them
pdf('Nichols_correlation_simple.pdf')
plot(clean_Sub_t_New_Data[,i],clean_Sub_t_New_Data[,j])











