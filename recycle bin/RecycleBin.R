#This is the Recycle Bin


#2nd kind of data: more than 1 number: (the smallest being 2 numbers. Eg. x<-matrix(c(1,2)))
#Might not need this part. Everything can be dealt using the 1st kind

#genotypeY+condY <-(-1.51)
##Should add a row and a column

#genotype + existing cond <-7.8
##Should add a row =>rbind()
sec_sec<-rbind(c(1.1,3.4),c(5.5,1.4))
colnames(sec_sec)<-c("Cond1","Cond2")
rownames(sec_sec)<-c("ECK2-2-1","ECK2-2-2")

rbind(small,sec_sec)

#exising genotype + condZ <-(-8.7)
##Should add a column => cbind()
sec_3rd<-rbind(c(1.1,3.4),c(5.5,1.4))
colnames(sec_3rd)<-c("Cond2-3-1","Cond2-3-2")
rownames(sec_3rd)<-c("ECK1","ECK2")

cbind(small,sec_3rd)


#existing genotype + existing cond
##Replicate => return an error and don't integrate
sec_4th<-rbind(c(1.1,3.4),c(5.5,1.4))
colnames(sec_4th)<-c("Cond1","Cond2")
rownames(sec_4th)<-c("ECK1","ECK2")


#Build up...
'
small<-rbind(c(1.1,3.4),c(5.5,1.4))
colnames(small)<-c("Cond1","Cond2")
rownames(small)<-c("ECK1","ECK2")
'

#Recycle Bin
'
SQL_matrix<-matrix(,0,3)
colnames(SQL_matrix)<-c(y_attribute,x_attribute,data_attribute)
rbind(SQL_matrix,t(as.matrix(c("ECK","UV",1))))
'



#Rotate function in dendextend can rotate according to a specified order (It tries to do so as much as possible)
#The ordering doesn't always work
par(mfrow = c(1,2))
All_Data[1:10,1:10] %>% dist() %>% hclust() %>% plot()
All_Data[1:10,1:10] %>% dist() %>% hclust() %>% dendextend::rotate(row.names(All_Data)[1:10]) %>% plot()


par(mfrow = c(1,2))
row<-50:100
col<-1:10
All_Data[row,col] %>% dist() %>% hclust() %>% dendextend::rotate(row.names(All_Data)[row]) %>% plot()
Binary_Data[row,col] %>% dist() %>% hclust() %>% dendextend::rotate(row.names(All_Data)[row]) %>% plot()


dendextend::cor_cophenetic(All_Data %>% dist() %>% hclust() %>% as.dendrogram(), Binary_Data %>% dist() %>% hclust() %>% as.dendrogram())
#Result: 0.6909157



