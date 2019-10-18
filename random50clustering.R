#Goal: get a heat map from randomly sampled 50*50 matrix from Nichols' data



#Load the heat map package
library(circlize)
library(ComplexHeatmap)



set.seed(1)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-All_Data[randomrows,randomcols]
pdf('random50clustering-1.pdf')
Heatmap(Sub_Data,col=colorRamp2(c(-4,0,4),c("blue","white","red")))
dev.off()


set.seed(2)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-All_Data[randomrows,randomcols]
pdf('random50clustering-2.pdf')
Heatmap(Sub_Data,col=colorRamp2(c(-4,0,4),c("blue","white","red")))
dev.off()

set.seed(3)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-All_Data[randomrows,randomcols]
pdf('random50clustering-3.pdf')
Heatmap(Sub_Data,col=colorRamp2(c(-4,0,4),c("blue","white","red")))
dev.off()

set.seed(4)
randomrows<-sample(1:3979,50)
randomcols<-sample(1:324,50)
Sub_Data<-All_Data[randomrows,randomcols]
pdf('random50clustering-4.pdf')
Heatmap(Sub_Data,col=colorRamp2(c(-4,0,4),c("blue","white","red")))
dev.off()