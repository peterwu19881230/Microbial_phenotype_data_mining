#Goal: Convert all data into binary (-1,0,1) using R


#This is not the most elegant way
matrix<-matrix(,3979,324)
rownames(matrix)<-rownames(All_Data)
colnames(matrix)<-colnames(All_Data)
matrix[All_Data>-3&All_Data<3]<-1 #Category 1: not significant
matrix[All_Data>=3]<-3 #Category 2: positive
matrix[All_Data<=-3]<-2  #Category 3: negative



write.csv(matrix, file = "categorical_Nichols.csv")


