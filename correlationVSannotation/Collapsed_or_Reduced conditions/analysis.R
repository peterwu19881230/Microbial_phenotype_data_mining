#This following graph gives PCC and different MHDs. It also tries to solve this problem: There are many same distances in MHD that should have been ranked differently
listOfResults=pwy_NAimputed.listOfResults


index_MHD3=listOfResults[[2]]$MHD3 %>% table %>% cumsum 


MHD3_NAnotImputed=strain1strain2.hammingNo0_corrected.samePWY_Annot$`Hamming Distance_no0_corrected`
reversed_table=MHD3_NAnotImputed %>% table %>% as.numeric #Have to use as.numeric to prevent furthur problems ("table" class doesn't do well with subsetting like [last_element:1])
index_MHD3_NAnotImputed=reversed_table[length(reversed_table):1] %>% cumsum #Had to reverse the order of this MHD because this is the earlier version where MHD doesn't really stand for distance
MHD3_NAnotImputed_cumsum=strain1strain2.hammingNo0_corrected.samePWY_Annot$bi.pwy.annot %>% cumsum

index_MHD4=listOfResults[[4]]$MHD4 %>% table %>% cumsum 


index_collapsed_MHD3=pwy_collapsedCond.MHD3[[1]]$MHD3 %>% table %>% cumsum
index_140cond_MHD3=nonRedun.140cond.Ternary_Data_MHD3[[1]]$MHD3 %>% table %>% cumsum 

index_collapsed_MHD4=pwy_collapsedCond.MHD4[[1]]$MHD4 %>% table %>% cumsum


samples=seq(from=0,to=7914231,by=1989); samples[1]=1 #change samples to 1:10000 to get the first 10000 results
#samples=1:dim(listOfResults[[1]])[1]
slopeForControl=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+geom_line(data=data.frame(no=samples,CumSum=listOfResults[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[3]]$cumsum[samples]),aes(no,CumSum,color="spearman"))+
  
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=listOfResults[[2]]$cumsum[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD3"))+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=listOfResults[[2]]$cumsum[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD3"))+
  
  #geom_point(data=data.frame(no=index_MHD3_NAnotImputed[index_MHD3_NAnotImputed<=max(samples)],CumSum=MHD3_NAnotImputed_cumsum[index_MHD3_NAnotImputed]),aes(no,CumSum,color="MHD_NAnotImputed"))+
  #geom_line(data=data.frame(no=index_MHD3_NAnotImputed[index_MHD3_NAnotImputed<=max(samples)],CumSum=MHD3_NAnotImputed_cumsum[index_MHD3_NAnotImputed]),aes(no,CumSum,color="MHD_NAnotImputed"))+
  
  #geom_point(data=data.frame(no=index_MHD4[index_MHD4<=max(samples)],CumSum=listOfResults[[4]]$cumsum[index_MHD4[index_MHD4<=max(samples)]]),aes(no,CumSum,color="MHD4"))+ 
  #geom_line(data=data.frame(no=index_MHD4[index_MHD4<=max(samples)],CumSum=listOfResults[[4]]$cumsum[index_MHD4[index_MHD4<=max(samples)]]),aes(no,CumSum,color="MHD4"))+
  #These 2 lines cause problems because sometimes index_MHD4[index_MHD4<=max(samples)] doesn't catch any index (eg. for MHD4=0 there are 176504 strain pairs)
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=pwy_collapsedCond.MHD3[[1]]$cumsum[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD3_collapsed"))+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=pwy_collapsedCond.MHD3[[1]]$cumsum[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD3_collapsed"))+
  
  geom_point(data=data.frame(no=index_140cond_MHD3[index_140cond_MHD3<=max(samples)],CumSum=nonRedun.140cond.Ternary_Data_MHD3[[1]]$cumsum[index_140cond_MHD3[index_140cond_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD3 140 less correlated condition"))+
  geom_line(data=data.frame(no=index_140cond_MHD3[index_140cond_MHD3<=max(samples)],CumSum=nonRedun.140cond.Ternary_Data_MHD3[[1]]$cumsum[index_140cond_MHD3[index_140cond_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD3 140 less correlated condition"))+
  
  #geom_point(data=data.frame(no=index_collapsed_MHD4[index_collapsed_MHD4<=max(samples)],CumSum=pwy_collapsedCond.MHD4[[1]]$cumsum[index_collapsed_MHD4[index_collapsed_MHD4<=max(samples)]]),aes(no,CumSum,color="MHD4_collapsed"))+
  #geom_line(data=data.frame(no=index_collapsed_MHD4[index_collapsed_MHD4<=max(samples)],CumSum=pwy_collapsedCond.MHD4[[1]]$cumsum[index_collapsed_MHD4[index_collapsed_MHD4<=max(samples)]]),aes(no,CumSum,color="MHD4_collapsed"))+
  #These 2 lines cause problems because sometimes index_MHD4[index_MHD4<=max(samples)] doesn't catch any index (eg. for MHD4=0 there are 176504 strain pairs)
  
  geom_abline(intercept=,slope=slopeForControl,linetype='dotted',color="black")+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="black")+ #positive control
  scale_y_continuous(expand = c(0,0),label=comma)+scale_x_continuous(expand = c(0,0),label=comma)+
  labs(x = "Low--Distance--High",y="No. of being in the same pwy",aesthetic='custom text')+
  scale_color_discrete(name="Distance")+
  theme_minimal()
