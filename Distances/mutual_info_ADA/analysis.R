load("Data/mutual_info_ternaryData_table.RData")
mutual_info_ternaryData_table=as.data.frame(mutual_info_ternaryData_table)
strain1strain2.samePWY_Annot_mutual_info=cbind(strain1strain2.samePWY_Annot,mutual_info=mutual_info_ternaryData_table[,3])
ordered_strain1strain2.samePWY_Annot_mutual_info=strain1strain2.samePWY_Annot_mutual_info[order(strain1strain2.samePWY_Annot_mutual_info$mutual_info,decreasing=T),]
mutual_info_cumsum=cumsum(ordered_strain1strain2.samePWY_Annot_mutual_info$bi.pwy.annot)

listOfResults=pwy_NAimputed.listOfResults

index_MHD3=listOfResults[[2]]$MHD3 %>% table %>% cumsum 

index_collapsed_MHD3=pwy_collapsedCond.MHD3[[1]]$MHD3 %>% table %>% cumsum 


#samples=seq(from=0,to=7914231,by=1989); samples[1]=1 #change samples to 1:10000 to get the first 10000 results
samples=1:4000
slopeForControl=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC based distance"))+
  geom_line(data=data.frame(no=samples,CumSum=mutual_info_cumsum[samples]),aes(no,CumSum,color="Mutual Info based distance"))+
  
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=listOfResults[[2]]$cumsum[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"))+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=listOfResults[[2]]$cumsum[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"))+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=pwy_collapsedCond.MHD3[[1]]$cumsum[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"))+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=pwy_collapsedCond.MHD3[[1]]$cumsum[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"))+
  
  geom_abline(intercept=,slope=slopeForControl,linetype='dotted',color="#D76E6E",size=1)+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="#1F89C2",size=1)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of pairs that have the same pathway annotation(s)",aesthetic='custom text')+
  scale_color_discrete(name="Distance")
