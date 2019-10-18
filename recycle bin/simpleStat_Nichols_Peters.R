#Goal: Do unit testing of correlation coefficient using part of Nichols' data



#Import the data from .csv
Nichols<-read.csv("allData_rows_fixed.csv",row.names=1,header = TRUE)
Peters<-read.csv("PetersAllKnockout.csv",row.names=1,header = TRUE)



#Convert them to double so mean() and sd() will work
Nichols<-data.matrix(Nichols)
Peters<-data.matrix(Peters)


#mean and sd of Nichols' data
mean(Nichols, na.rm=TRUE)
sd(Nichols, na.rm=TRUE)

#mean and sd of Peters' data
mean(Peters,na.rm=TRUE)
sd(Peters,na.rm=TRUE)



#Output the histogram to pdf files
pdf('Nichols.pdf')
hist(Nichols,xlab ='Normalized growth scores')

pdf('Peters.pdf')
hist(Peters,xlab ='Normalized growth scores')















