# CorrelationVSAnnotation analysis using KEGG pathways
# 1. For everything (3979 strains used in Nichols')
# 2. For Nichols' strains that have been annotated with KEGG modules



# 1. For everything (3979 strains used in Nichols')

modifiedSameAnnot=numeric(dim(MHD3distance_table)[1]) #An all-zero vector
for(mhd in unique(MHD3distance_table$Distance)) {
  index=(MHD3distance_table$Distance==mhd)
  modifiedSameAnnot[index]=mean(MHD3distance_table$sameORnot[index])
}
modifiedCumsum=cumsum(modifiedSameAnnot)

listOfResults=kegg_quantitative
slopeForControl=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
samples=seq(from=0,to=7914231,by=1989); samples[1]=1
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+
  #geom_line(data=data.frame(no=samples,CumSum=cumSum_MHDdistance_table_TF[samples]),aes(no,CumSum,color="MHD"))+
  #geom_line(data=data.frame(no=samples,CumSum=cumSum_MHD2distance_table_TF[samples]),aes(no,CumSum,color="MHD2"))+
  geom_line(data=data.frame(no=samples,CumSum=cumSum_MHD3distance_table_TF[samples]),aes(no,CumSum,color="MHD3"))+
  geom_line(data=data.frame(no=samples,CumSum=modifiedCumsum[samples]),aes(no,CumSum,color="MHD3_modifiedCumsum"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[2]]$cumsum[samples]),aes(no,CumSum,color="PCC_squared"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[3]]$cumsum[samples]),aes(no,CumSum,color="euclidean"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[4]]$cumsum[samples]),aes(no,CumSum,color="maximum"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[5]]$cumsum[samples]),aes(no,CumSum,color="manhattan"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[6]]$cumsum[samples]),aes(no,CumSum,color="canberra"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[7]]$cumsum[samples]),aes(no,CumSum,color="binary"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[8]]$cumsum[samples]),aes(no,CumSum,color="minkowski"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[9]]$cumsum[samples]),aes(no,CumSum,color="spearman"))+  
  geom_abline(intercept=,slope=slopeForControl,linetype='dotted',color="black")+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="black")+ #positive control
  scale_y_continuous(expand = c(0,0),label=comma)+scale_x_continuous(expand = c(0,0),label=comma)+
  labs(x = "Low--Distance--High",y="No. of being in the same KEGG modules",aesthetic='custom text')+
  scale_color_discrete(name="Distance")







# 2. For Nichols' strains that have been annotated with KEGG modules

modifiedSameAnnot=numeric(dim(MHD3distanceInModule_table)[1]) #An all-zero vector
for(mhd in unique(MHD3distanceInModule_table$Distance)) {
  index=(MHD3distanceInModule_table$Distance==mhd)
  modifiedSameAnnot[index]=mean(MHD3distanceInModule_table$sameORnot[index])
}
modifiedCumsum=cumsum(modifiedSameAnnot)

listOfResults=keggInModules_quantitative
slopeForControl=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
samples=1:dim(listOfResults[[1]])[1]
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+
  #geom_line(data=data.frame(no=samples,CumSum=cumSum_MHDdistanceInModule_table_TF[samples]),aes(no,CumSum,color="MHD"))+
  #geom_line(data=data.frame(no=samples,CumSum=cumSum_MHD2distanceInModule_table_TF[samples]),aes(no,CumSum,color="MHD2"))+
  geom_line(data=data.frame(no=samples,CumSum=cumSum_MHD3distanceInModule_table_TF[samples]),aes(no,CumSum,color="MHD3"))+
  geom_line(data=data.frame(no=samples,CumSum=modifiedCumsum[samples]),aes(no,CumSum,color="MHD3_modifiedCumsum"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[2]]$cumsum[samples]),aes(no,CumSum,color="PCC_squared"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[3]]$cumsum[samples]),aes(no,CumSum,color="euclidean"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[4]]$cumsum[samples]),aes(no,CumSum,color="maximum"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[5]]$cumsum[samples]),aes(no,CumSum,color="manhattan"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[6]]$cumsum[samples]),aes(no,CumSum,color="canberra"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[7]]$cumsum[samples]),aes(no,CumSum,color="binary"))+
  #geom_line(data=data.frame(no=samples,CumSum=listOfResults[[8]]$cumsum[samples]),aes(no,CumSum,color="minkowski"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[9]]$cumsum[samples]),aes(no,CumSum,color="spearman"))+  
  geom_abline(intercept=,slope=slopeForControl,linetype='dotted',color="black")+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="black")+ #positive control
  scale_y_continuous(expand = c(0,0),label=comma)+scale_x_continuous(expand = c(0,0),label=comma)+
  labs(x = "Low--Distance--High",y="No. of being in the same KEGG modules",aesthetic='custom text')+
  scale_color_discrete(name="Distance")






