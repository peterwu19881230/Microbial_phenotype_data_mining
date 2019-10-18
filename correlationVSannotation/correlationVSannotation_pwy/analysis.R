##Function to do lineplots of all & the first 10000 & the first 100 (Would take too long to plot the whole pcc)
compareDist_cummulativeFraction=function(quantitative_TF,quantitative_coef,qualitative_TF,qualitative_coef){
  par(mfrow=c(1,2))
  
  fraction=cumsum(quantitative_TF)/(1:length(quantitative_TF))
  fraction2=cumsum(qualitative_TF)/(1:length(qualitative_TF))
  
  ##This is to do all:
  #plot(main="All coef V.S. pwy annotation in Nichols' (pseudo-cummulative)",1:length(fraction),fraction,type="l",xlab="Order of the decreasing coef",ylab="No. of the same pwy annotation / No. of coef until the  cutoff coef",ylim=c(0,1))
  #lines(1:length(fraction2),fraction2,type="l",col="red",ylim=c(0,1))
  
  plot(main="First 10000 coef V.S. pwy annotation in Nichols' (pseudo-cummulative)",1:10000,fraction[1:10000],type="l",xlab="Order of the decreasing coef",ylab="No. of the same pwy annotation / No. of coef until the  cutoff coef",ylim=c(0,1)) 
  
  text(c(1,2000,4000,6000,8000,10000),
       fraction[c(1,2000,4000,6000,8000,10000)],pos=3,
       labels=round(c(quantitative_coef[1],quantitative_coef[2000],quantitative_coef[4000],quantitative_coef[6000],quantitative_coef[8000],quantitative_coef[10000]),4),
       cex=0.5) #"3" in pos=3 means "above" 
  lines(1:10000,fraction2[1:10000],type="l",col="red",ylim=c(0,1))
  text(c(1,2000,4000,6000,8000,10000),
       fraction2[c(1,2000,4000,6000,8000,10000)],pos=3,
       labels=round(c(qualitative_coef[1],qualitative_coef[2000],qualitative_coef[4000],qualitative_coef[6000],qualitative_coef[8000],qualitative_coef[10000]),4),cex=0.5,col="red") #"3" in pos=3 means "above"
  
  
  plot(main="Fig. 1(c): First 100 quantitative_coef V.S. pwy annotation in Nichols' (pseudo-cummulative)",1:100,fraction[1:100],type="l",xlab="Order of the decreasing coef",ylab="No. of the same pwy annotation / No. of coef until the  cutoff coef",ylim=c(0,1))
  text(c(1,20,40,60,80,100),
       fraction[c(1,20,40,60,80,100)],pos=3,
       labels=round(c(quantitative_coef[1],quantitative_coef[20],quantitative_coef[40],quantitative_coef[60],quantitative_coef[80],quantitative_coef[100]),4),
       cex=0.5) #"3" in pos=3 means "above"
  lines(1:100,fraction2[1:100],type="l",col="red",ylim=c(0,1))
  text(c(1,20,40,60,80,100),
       fraction2[c(1,20,40,60,80,100)],pos=3,
       labels=round(c(qualitative_coef[1],qualitative_coef[20],qualitative_coef[40],qualitative_coef[60],qualitative_coef[80],qualitative_coef[100]),4),
       cex=0.5,col="red") #"3" in pos=3 means "above" 
  
  
  par(mfrow=c(1,1))
}


#All gene pairs using Nichols' fig.S1 pwy
compareDist_cummulativeFraction(pcc_TF,pcc_coef,MHDnoPunishment_TF,MHDnoPunishment_coef)
compareDist_cummulativeFraction(pcc_TF,pcc_coef,MHD_TF,MHD_coef)
compareDist_cummulativeFraction(inPwy.pcc_TF,inPwy.pcc_coef,inPwy.MHDnoPunishment_TF,inPwy.MHDnoPunishment_coef)
compareDist_cummulativeFraction(inPwy.pcc_TF,inPwy.pcc_coef,inPwy.MHD_TF,inPwy.MHD_coef)


#This following graph gives PCC and different MHDs. It also tries to solve this problem: Many same distances in MHD that shouldn't have been ranked differently
listOfResults=pwy_NAimputed.listOfResults

modifiedSameAnnot=numeric(dim(listOfResults[[2]])[1]) #An all-zero vector
for(mhd in unique(listOfResults[[2]]$MHD3)) {
  index=(listOfResults[[2]]$MHD3==mhd)
  modifiedSameAnnot[index]=mean(listOfResults[[2]]$sameORnot[index])
}
modifiedCumsum=cumsum(modifiedSameAnnot)


modifiedSameAnnot_NAnotImputed=numeric(dim(strain1strain2.hammingNo0_corrected.samePWY_Annot)[1]) #An all-zero vector
for(mhd in unique(strain1strain2.hammingNo0_corrected.samePWY_Annot$`Hamming Distance_no0_corrected`)) {
  index=(strain1strain2.hammingNo0_corrected.samePWY_Annot$`Hamming Distance_no0_corrected`==mhd)
  modifiedSameAnnot_NAnotImputed[index]=mean(strain1strain2.hammingNo0_corrected.samePWY_Annot$bi.pwy.annot[index])
}
modifiedCumsum_NAnotImputed=cumsum(modifiedSameAnnot_NAnotImputed)


samples=seq(from=0,to=7914231,by=1989); samples[1]=1 #change samples to 1:10000 to get the first 10000 results
#samples=1:dim(listOfResults[[1]])[1]
slopeForControl=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+geom_line(data=data.frame(no=samples,CumSum=listOfResults[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[3]]$cumsum[samples]),aes(no,CumSum,color="spearman"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[2]]$cumsum[samples]),aes(no,CumSum,color="MHD"))+
  geom_line(data=data.frame(no=samples,CumSum=modifiedCumsum[samples]),aes(no,CumSum,color="MHD_modifiedCumsum"))+
  geom_line(data=data.frame(no=samples,CumSum=cumsum(MHD_TF)[samples]),aes(no,CumSum,color="MHD_NAnotImputed"))+
  geom_line(data=data.frame(no=samples,CumSum=modifiedCumsum_NAnotImputed[samples]),aes(no,CumSum,color="MHD_NAnotImputed_modifiedCumsum"))+
  geom_abline(intercept=,slope=slopeForControl,linetype='dotted',color="black")+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="black")+ #positive control
  scale_y_continuous(expand = c(0,0),label=comma)+scale_x_continuous(expand = c(0,0),label=comma)+
  labs(x = "Low--Distance--High",y="No. of being in the same pwy",aesthetic='custom text')+
  scale_color_discrete(name="Distance")


#(!)Make sure this barplot is what we want then finish the remaining code (Only genes involved in Nichols' fig.S1 pwy, ptcomplex, KEGG,...etc)
#barplot
uniqueMHD=unique(listOfResults[[2]]$MHD3) %>% length
index=listOfResults[[2]]$MHD3 %>% table %>% cumsum 
PCC_bars=listOfResults[[1]]$cumsum[index]
MHD_bars=listOfResults[[2]]$cumsum[index]
MHD4_bars=listOfResults[[4]]$cumsum[index]
spearman_bars=listOfResults[[3]]$cumsum[index]
negative_control=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
positive_control=ifelse(index<=sum(listOfResults[[1]]$sameORnot),index,sum(listOfResults[[1]]$sameORnot))
  

subset=1:71 #maximum: 1:97. 71 is the last one before the negative control starts climbing


##side-by-side (didn't use the spearman. Otherwise it will be too messy)
#dat=data.frame(MHD_distance=rep(unique(listOfResults[[2]]$MHD3)[subset],2) %>% as.character,PCC_MHD=c(PCC_bars[subset],MHD_bars[subset]),name=c(rep("PCC",uniqueMHD)[subset],rep("MHD",uniqueMHD)[subset]))
#ggplot(data=dat,aes(x=MHD_distance,y=PCC_MHD,fill=name))+
#  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal()+scale_fill_brewer(palette="Blues")+
#  geom_line(data=data.frame(x=dat$MHD_distance[subset] %>% as.character,y=negative_control*index[subset],name="negative"),aes(x=x,y=y,group="negative"),color="red")+
#  geom_line(data=data.frame(x=dat$MHD_distance[subset] %>% as.character,y=positive_control[subset],name="positive"),aes(x=x,y=y,group="positive"),color="blue")

##PCC is plotted as a line
##Ref for text size&orientation adjustment: https://stackoverflow.com/questions/13297995/changing-font-size-and-direction-of-axes-text-in-ggplot2
dat=data.frame(MHD_distance=unique(listOfResults[[2]]$MHD3)[subset] %>% as.character,MHD=MHD_bars[subset],name=rep("MHD",uniqueMHD)[subset])
ggplot(data=dat,aes(x=MHD_distance,y=MHD,fill=name))+
  geom_bar(stat="identity", color="black")+theme_minimal()+theme(axis.text.y = element_text(size=20),axis.text.x = element_text(size=14))+scale_fill_brewer(palette="Blues")+
  geom_line(data=data.frame(x=dat$MHD_distance[subset] %>% as.character,y=PCC_bars[subset],name="PCC"),aes(x=x,y=y,group="PCC"),size=1,color="purple")+
  geom_line(data=data.frame(x=dat$MHD_distance[subset] %>% as.character,y=spearman_bars[subset],name="spearman"),aes(x=x,y=y,group="spearman"),size=1,color="pink")+
  geom_line(data=data.frame(x=dat$MHD_distance[subset] %>% as.character,y=negative_control*index[subset],name="negative"),aes(x=x,y=y,group="negative"),linetype="dashed",size=1,color="red",alpha=0.7)+
  geom_line(data=data.frame(x=dat$MHD_distance[subset] %>% as.character,y=positive_control[subset],name="positive"),aes(x=x,y=y,group="positive"),linetype="dashed",size=1,color="#225BC6",alpha=0.7)


##Using MHD4. PCC is plotted as a line
uniqueMHD4=unique(listOfResults[[4]]$MHD4) %>% length

###==== table() doesn't seem to work for NaN. Eg. table(c("NaN",1,2,3,4)) doesn't give the count for NaN. I have to change NaN to another value 
listOfResults[[4]]$MHD4=listOfResults[[4]]$MHD4 %>% as.character
listOfResults[[4]]$MHD4[listOfResults[[4]]$MHD4=="NaN"]="nan"
###====

index=listOfResults[[4]]$MHD4 %>% table %>% cumsum 
PCC_bars=listOfResults[[1]]$cumsum[index]
MHD4_bars=listOfResults[[4]]$cumsum[index]
spearman_bars=listOfResults[[3]]$cumsum[index]
negative_control=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
positive_control=ifelse(index<=sum(listOfResults[[1]]$sameORnot),index,sum(listOfResults[[1]]$sameORnot))

subset=1:(listOfResults[[4]]$MHD4 %>% unique %>% length) #maximum: 195 (194 distances + NA =195) 

##Ref for text size&orientation adjustment: https://stackoverflow.com/questions/13297995/changing-font-size-and-direction-of-axes-text-in-ggplot2
dat=data.frame(MHD4_distance=unique(listOfResults[[4]]$MHD4)[subset] %>% as.character,MHD4=MHD4_bars[subset],name=rep("MHD4",uniqueMHD)[subset])
ggplot(data=dat,aes(x=MHD4_distance,y=MHD4,fill=name))+
  geom_bar(stat="identity", color="black")+theme_minimal()+theme(axis.text.y = element_text(size=20),axis.text.x = element_text(size=14))+scale_fill_brewer(palette="Blues")+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=PCC_bars[subset],name="PCC"),aes(x=x,y=y,group="PCC"),size=1,color="purple")+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=spearman_bars[subset],name="spearman"),aes(x=x,y=y,group="spearman"),size=1,color="pink")+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=negative_control*index[subset],name="negative control"),aes(x=x,y=y,group="negative control"),linetype="dashed",size=1,color="red",alpha=0.7)+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=positive_control[subset],name="positive control"),aes(x=x,y=y,group="positive control"),linetype="dashed",size=1,color="#225BC6",alpha=0.7)


##Line plot
ggplot()+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=PCC_bars[subset],name="PCC"),aes(x=x,y=y,group="PCC"),size=1,color="purple")+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=spearman_bars[subset],name="spearman"),aes(x=x,y=y,group="spearman"),size=1,color="pink",alpha=0.7)+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=negative_control*index[subset],name="negative control"),aes(x=x,y=y,group="negative control"),linetype="dashed",size=1,color="red",alpha=0.7)+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=positive_control[subset],name="positive control"),aes(x=x,y=y,group="positive control"),linetype="dashed",size=1,color="#225BC6",alpha=0.7)+
  geom_line(data=data.frame(x=dat$MHD4_distance[subset] %>% as.character,y=MHD4_bars[subset],name="MHD4"),aes(x=x,y=y,group="MHD4"),size=1,color="black")


##conclusion for MHD4: 
## When distance=0 (smallest distance. Total No.=1211), I see that there are lots of them that don't have the same annotation(s)



#This following graph gives PCC and different MHDs. It also tries to solve this problem: There are many same distances in MHD that should have been ranked differently
MHD3_NAnotImputed=strain1strain2.hammingNo0_corrected.samePWY_Annot$`Hamming Distance_no0_corrected`
reversed_table=MHD3_NAnotImputed %>% table %>% as.numeric #Have to use as.numeric to prevent furthur problems ("table" class doesn't do well with subsetting like [last_element:1])
index_MHD3_NAnotImputed=reversed_table[length(reversed_table):1] %>% cumsum #Had to reverse the order of this MHD because this is the earlier version where MHD doesn't really stand for distance
MHD3_NAnotImputed_cumsum=strain1strain2.hammingNo0_corrected.samePWY_Annot$bi.pwy.annot %>% cumsum

index_MHD3=listOfResults[[2]]$MHD3 %>% table %>% cumsum
index_MHD4=listOfResults[[4]]$MHD4 %>% table %>% cumsum
index_MHD5=listOfResults[[5]]$MHD5 %>% table %>% cumsum 


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
  
  geom_point(data=data.frame(no=index_MHD4[index_MHD4<=max(samples)],CumSum=listOfResults[[4]]$cumsum[index_MHD4[index_MHD4<=max(samples)]]),aes(no,CumSum,color="MHD4"))+ 
  geom_line(data=data.frame(no=index_MHD4[index_MHD4<=max(samples)],CumSum=listOfResults[[4]]$cumsum[index_MHD4[index_MHD4<=max(samples)]]),aes(no,CumSum,color="MHD4"))+
  #These 2 lines cause problems because sometimes index_MHD4[index_MHD4<=max(samples)] doesn't catch any index (eg. for MHD4=0 there are 176504 strain pairs)
  
  geom_point(data=data.frame(no=index_MHD5[index_MHD5<=max(samples)],CumSum=listOfResults[[5]]$cumsum[index_MHD5[index_MHD5<=max(samples)]]),aes(no,CumSum,color="MHD5"))+ 
  geom_line(data=data.frame(no=index_MHD5[index_MHD5<=max(samples)],CumSum=listOfResults[[5]]$cumsum[index_MHD5[index_MHD5<=max(samples)]]),aes(no,CumSum,color="MHD5"))+
  #These 2 lines cause problems because sometimes index_MHD5[index_MHD5<=max(samples)] doesn't catch any index (eg. for MHD4=0 there are 176504 strain pairs)
  
  geom_abline(intercept=,slope=slopeForControl,linetype='dotted',color="black")+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="black")+ #positive control
  scale_y_continuous(expand = c(0,0),label=comma)+scale_x_continuous(expand = c(0,0),label=comma)+
  labs(x = "Low--Distance--High",y="No. of being in the same pwy",aesthetic='custom text')+
  scale_color_discrete(name="Distance")


#Conclusioin: 
##1. From the first 200,000 I know that MHD4 doesn't work well  
##Note: from this I know that in all strain pairs where MHD4=0, only a tiny portion of them have the same annotations: 
##((listOfResults[[4]]$MHD4==0 & listOfResults[[4]]$sameORnot==1) %>% sum(na.rm=T) ) / ((listOfResults[[4]]$MHD4==0) %>% sum(na.rm=T)) # 0.007620224
##2. From the first 10,000 I know that MHD5 doesn't work well
##Note: When distance=0 (smallest distance. Total No.=1345), I see that there are lots of them that don't have the same annotation(s)
##Note: from this I know that in all strain pairs where MHD5=0, only a tiny portion of them have the same annotations: 
##((listOfResults[[5]]$MHD5==0 & listOfResults[[5]]$sameORnot==1) %>% sum(na.rm=T) ) / ((listOfResults[[5]]$MHD5==0) %>% sum(na.rm=T)) # 0.03385632












#Only genes involved in Nichols' fig.S1 pwy

listOfResults=inpwy_NAimputed.listOfResults

modifiedSameAnnot=numeric(dim(listOfResults[[2]])[1]) #An all-zero vector
for(mhd in unique(listOfResults[[2]]$MHD3)) {
  index=(listOfResults[[2]]$MHD3==mhd)
  modifiedSameAnnot[index]=mean(listOfResults[[2]]$sameORnot[index])
}
modifiedCumsum=cumsum(modifiedSameAnnot)

strainsInPwy=unique(c(pwy.pcc$strain1,pwy.pcc$strain2))
table=strain1strain2.hammingNo0_corrected.samePWY_Annot[strain1strain2.hammingNo0_corrected.samePWY_Annot$strain1 %in% strainsInPwy &
                                                        strain1strain2.hammingNo0_corrected.samePWY_Annot$strain2 %in% strainsInPwy,]
modifiedSameAnnot_NAnotImputed=numeric(dim()[1]) #An all-zero vector
for(mhd in unique(table$`Hamming Distance_no0_corrected`)) {
  index=(table$`Hamming Distance_no0_corrected`==mhd)
  modifiedSameAnnot_NAnotImputed[index]=mean(table$bi.pwy.annot[index])
}
modifiedCumsum_NAnotImputed=cumsum(modifiedSameAnnot_NAnotImputed)

#samples=seq(from=0,to=7914231,by=1989); samples[1]=1 #change samples to 1:10000 to get the first 10000 results
samples=1:dim(listOfResults[[1]])[1]
slopeForControl=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+geom_line(data=data.frame(no=samples,CumSum=listOfResults[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[2]]$cumsum[samples]),aes(no,CumSum,color="MHD"))+
  geom_line(data=data.frame(no=samples,CumSum=modifiedCumsum[samples]),aes(no,CumSum,color="MHD_modifiedCumsum"))+
  geom_line(data=data.frame(no=samples,CumSum=cumsum(inPwy.MHD_TF)[samples]),aes(no,CumSum,color="MHD_NAnotImputed"))+
  geom_line(data=data.frame(no=samples,CumSum=modifiedCumsum_NAnotImputed[samples]),aes(no,CumSum,color="MHD_NAnotImputed_modifiedCumsum"))+
  geom_abline(intercept=,slope=slopeForControl,linetype='dotted',color="black")+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="black")+ #positive control
  scale_y_continuous(expand = c(0,0),label=comma)+scale_x_continuous(expand = c(0,0),label=comma)+
  labs(x = "Low--Distance--High",y="No. of being in the same pwy",aesthetic='custom text')+
  scale_color_discrete(name="Distance")

