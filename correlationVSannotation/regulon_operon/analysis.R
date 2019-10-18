#For regulons

index_MHD3=regulonAnnotation[[3]]$distance %>% table %>% cumsum 

index_collapsed_MHD3=regulonAnnotation[[4]]$distance %>% table %>% cumsum

samples=1:4000
slopeForControl=sum(regulonAnnotation[[1]]$sameORnot)/length(regulonAnnotation[[1]]$sameORnot)
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,CumSum=regulonAnnotation[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC based distance"))+
  geom_line(data=data.frame(no=samples,CumSum=regulonAnnotation[[2]]$cumsum[samples]),aes(no,CumSum,color="Spearman based distance"))+
  
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=regulonAnnotation[[3]]$cumsum[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"))+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=regulonAnnotation[[3]]$cumsum[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"))+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=regulonAnnotation[[4]]$cumsum[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"))+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=regulonAnnotation[[4]]$cumsum[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"))+
  
  geom_abline(intercept=0,slope=slopeForControl,linetype='dotted',color="#D76E6E",size=1)+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="#1F89C2",size=1)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of pairs that have the same pathway annotation(s)",aesthetic='custom text')+
  scale_color_discrete(name="Distance")








#For operons
index_MHD3=operonAnnotation[[3]]$distance %>% table %>% cumsum 

index_collapsed_MHD3=operonAnnotation[[4]]$distance %>% table %>% cumsum

samples=1:4000
slopeForControl=sum(operonAnnotation[[1]]$sameORnot)/length(operonAnnotation[[1]]$sameORnot)
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,CumSum=operonAnnotation[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC based distance"))+
  geom_line(data=data.frame(no=samples,CumSum=operonAnnotation[[2]]$cumsum[samples]),aes(no,CumSum,color="Spearman based distance"))+
  
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=operonAnnotation[[3]]$cumsum[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"))+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=operonAnnotation[[3]]$cumsum[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"))+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=operonAnnotation[[4]]$cumsum[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"))+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=operonAnnotation[[4]]$cumsum[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"))+
  
  geom_abline(intercept=0,slope=slopeForControl,linetype='dotted',color="#D76E6E",size=1)+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="#1F89C2",size=1)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of pairs that have the same pathway annotation(s)",aesthetic='custom text')+
  scale_color_discrete(name="Distance")

