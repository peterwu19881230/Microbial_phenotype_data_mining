#Goal: get a heat map from randomly sampled 50*50 matrix from Nichols' data


#Load the heat map package
library(circlize)
library(ComplexHeatmap)



#Binary
matrix<-Binary_Data
pdf('binary_random50clustering.pdf')

set.seed(1)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-matrix[randomrows,randomcols]
Heatmap(Sub_Data,name="1",col=colorRamp2(c(-1,0,1),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE,cluster_rows=FALSE,cluster_columns=FALSE)
Heatmap(Sub_Data,name="1_clustered",col=colorRamp2(c(-1,0,1),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE)


set.seed(2)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-matrix[randomrows,randomcols]
Heatmap(Sub_Data,name="2",col=colorRamp2(c(-1,0,1),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE,cluster_rows=FALSE,cluster_columns=FALSE)
Heatmap(Sub_Data,name="2_clustered",col=colorRamp2(c(-1,0,1),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE)


set.seed(3)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-matrix[randomrows,randomcols]
Heatmap(Sub_Data,name="3",col=colorRamp2(c(-1,0,1),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE,cluster_rows=FALSE,cluster_columns=FALSE)
Heatmap(Sub_Data,name="3_clustered",col=colorRamp2(c(-1,0,1),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE)


#set.seed(4)
#randomrows<-sample(1:3979,50)
#randomcols<-sample(1:324,50)
#Sub_Data<-matrix[randomrows,randomcols]
#Heatmap(Sub_Data,name="4",col=colorRamp2(c(-1,0,1),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE,cluster_rows=FALSE,cluster_columns=FALSE)
#Heatmap(Sub_Data,name="4_clustered",col=colorRamp2(c(-1,0,1),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE)

dev.off()


#Original

matrix<-All_Data
pdf('original_random50clustering.pdf')

set.seed(1)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-matrix[randomrows,randomcols]
Heatmap(Sub_Data,name="1",col=colorRamp2(c(-3,0,3),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE,cluster_rows=FALSE,cluster_columns=FALSE)
Heatmap(Sub_Data,name="1_clustered",col=colorRamp2(c(-3,0,3),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE)


set.seed(2)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-matrix[randomrows,randomcols]
Heatmap(Sub_Data,name="2",col=colorRamp2(c(-3,0,3),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE,cluster_rows=FALSE,cluster_columns=FALSE)
Heatmap(Sub_Data,name="2_clustered",col=colorRamp2(c(-3,0,3),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE)


set.seed(3)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-matrix[randomrows,randomcols]
Heatmap(Sub_Data,name="3",col=colorRamp2(c(-3,0,3),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE,cluster_rows=FALSE,cluster_columns=FALSE)
Heatmap(Sub_Data,name="3_clustered",col=colorRamp2(c(-3,0,3),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE)


#set.seed(4)
#randomrows<-sample(1:3979,50)
#randomcols<-sample(1:324,50)
#Sub_Data<-matrix[randomrows,randomcols]
#Heatmap(Sub_Data,name="4",col=colorRamp2(c(-3,0,3),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE,cluster_rows=FALSE,cluster_columns=FALSE)
#Heatmap(Sub_Data,name="4_clustered",col=colorRamp2(c(-3,0,3),c("blue","white","red")),show_row_names=FALSE,show_column_names = FALSE)

dev.off()
