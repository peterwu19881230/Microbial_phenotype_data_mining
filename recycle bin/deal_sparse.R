#Goal: test how sparse data can be integrated into the original matrix
#I might not need this...I can manually put multiple spread sheets into 1 big csv


#These are matrix for testing the function:

#create original matrix
small<-rbind(c(1.1,3.4),c(5.5,1.4))
colnames(small)<-c("Cond1","Cond2")
rownames(small)<-c("ECK1","ECK2")

#Another matrix (test if the function can add elemetns with new row and col names)
add<-rbind(c(-2,4.38),c(1.11,4.43))
colnames(add)<-c("Cond3","Cond4")
rownames(add)<-c("ECK3","ECK4")

#3rd matrix (test if the function can complain if there are duplicated values)
add2<-rbind(c(-1.9,4.33),c(1.2,4.5))
colnames(add2)<-c("Cond1","Cond2")
rownames(add2)<-c("ECK1","ECK2")




#Test if I can iterate several times:
mat3<-All_Data[1:10,1:10]
mat4<-All_Data[5:15,5:15]
mat5<-All_Data[20:30,20:30]
mat6<-All_Data[1000,300:301]
mat7<-All_Data[2500,323:324]


#single value shouldn't work. My function input should be a matrix



#Test overlaying matrices by subsetting All_Data. Combining these 2 takes about 20 sec on my mac
mat1<-All_Data[1:100,1:100]
mat2<-All_Data[90:200,90:200]






#I moved the function to function.R. It's called "matADDdf"
