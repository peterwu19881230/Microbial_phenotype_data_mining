#Correlation VS. annotation from all annotation sources and by all distances

sigPheno_strain1strain2_allAnnotations_allDistances=strain1strain2_allAnnotations_allDistances[strain1strain2_allAnnotations_allDistances$sigPhenoInBoth,]
dim(sigPheno_strain1strain2_allAnnotations_allDistances)

anyCoAnnotation=ifelse(rowSums(sigPheno_strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>0,
                       T,F)

cumsums=list()
distances=c("pcc","spearman","euclidean","euclidean_qualitative","euclidean_collapsedCond","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond")
for(distance in distances){
  cumsums[[distance]]=cumsum(anyCoAnnotation[order(sigPheno_strain1strain2_allAnnotations_allDistances[,distance])])
}

str(cumsums)


#plot it
index_MHD3=sigPheno_strain1strain2_allAnnotations_allDistances$mhd3 %>% table %>% cumsum 
index_collapsed_MHD3=sigPheno_strain1strain2_allAnnotations_allDistances$mhd3_collapsedCond %>% table %>% cumsum
slopeForControl=sum(anyCoAnnotation)/length(anyCoAnnotation)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work


#samples=seq(from=0,to=2440945,by=799); samples[1]=1
samples=1:4000
size=1
alpha=0.5

ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,CumSum=cumsums$pcc[samples]),aes(no,CumSum,color="PCC based distance"),size=size,alpha=alpha,linetype=2)+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$spearman[samples]),aes(no,CumSum,color="Spearman based distance"),size=size,alpha=alpha,linetype=3)+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$euclidean[samples]),aes(no,CumSum,color="Euclidean based distance"),size=size,alpha=alpha,linetype=4)+
  
  geom_line(data=data.frame(no=samples,CumSum=cumsums$euclidean_qualitative[samples]),aes(no,CumSum,color="Euclidean qualitative data"),size=size,alpha=alpha,linetype=5)+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$euclidean_collapsedCond[samples]),aes(no,CumSum,color="Euclidean collapsed condition"),size=size,alpha=alpha,linetype=6)+
  
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"),size=size*2,alpha=alpha)+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"),size=size,alpha=alpha,linetype=8)+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"),size=size*2,alpha=alpha)+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"),size=size,alpha=alpha,linetype=9)+
  
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi[samples]),aes(no,CumSum,color="MI based distance"),size=size,alpha=alpha,linetype=10)+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi_collapsedCond[samples]),aes(no,CumSum,color="MI based distance - collapsed condition"),size=size,alpha=alpha,linetype=11)+
  
  geom_abline(intercept=0,slope=slopeForControl,color="#D76E6E",size=1,alpha=alpha)+ #negative control
  geom_abline(intercept=0, slope=1,color="#1F89C2",size=1,alpha=alpha)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of pairs that have the same pathway annotation(s)",aesthetic='custom text')+
  scale_color_discrete(name="Distance")


#barplot at particular correlation cutoffs (for: "pcc","spearman","euclidean"): ranking of distance=c(100,200,500,1000,2000,3000)

##Check the corresponding pcc values
cutoffs=c(100,200,500,1000,2000,3000)
sort(1-sigPheno_strain1strain2_allAnnotations_allDistances$pcc,decreasing=T)[cutoffs]


p=list()
for(i in 1:length(cutoffs)){
  dat=data.frame(Distance=c("pcc","spearman","euclidean","random expectation"),value=c(cumsums$pcc[cutoffs[i]],cumsums$spearman[cutoffs[i]],cumsums$euclidean[cutoffs[i]],sum(anyCoAnnotation)/length(anyCoAnnotation)*cutoffs[i]))
  dat$Distance=factor(dat$Distance,levels=dat$Distance)#Do this to prevent automatic x axis reordering
  p[[i]]=ggplot(data=dat,aes(Distance,value))+theme_minimal()+geom_bar(stat="identity")+labs(y="Cummulative no. of co-annotations")+labs(x="")
}

gridExtra::grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],ncol=1)

dat=data.frame()
for(i in 1:length(cutoffs)){
  dat=rbind(dat,
            data.frame(Distance=c("pcc","spearman","euclidean","random expectation"),
                       value=c(cumsums$pcc[cutoffs[i]],cumsums$spearman[cutoffs[i]],cumsums$euclidean[cutoffs[i]],sum(anyCoAnnotation)/length(anyCoAnnotation)*cutoffs[i]),
                       cutoff=rep(i,4)))
}
dat$Distance=factor(dat$Distance,levels=c("pcc","spearman","euclidean","random expectation"))#Do this to prevent automatic x axis reordering

ggplot(data=dat,aes(Distance,value,fill=cutoff))+
  geom_bar(stat="identity")+
  labs(y="Cummulative no. of co-annotations")+
  scale_fill_gradient(labels=cutoffs,low="#5F06CD",high="#E0DFF7")
  

##The corresponding plot for "pcc","spearman","euclidean"

samples=1:4000
alpha=0.5
size=1

ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,PCC_based_distance=cumsums$pcc[samples]),aes(no,PCC_based_distance,color="PCC based distance"),size=size,alpha=alpha,linetype=1)+
  geom_line(data=data.frame(no=samples,Spearman_based_distance=cumsums$spearman[samples]),aes(no,Spearman_based_distance,color="Spearman based distance"),size=size,alpha=alpha,linetype=2)+
  geom_line(data=data.frame(no=samples,Euclidean_based_distance=cumsums$euclidean[samples]),aes(no,Euclidean_based_distance,color="Euclidean based distance"),size=size,alpha=alpha,linetype=3)+
  
  geom_abline(intercept=0,slope=slopeForControl,color="#D76E6E",size=1,alpha=alpha)+ #negative control
  geom_abline(intercept=0, slope=1,color="#1F89C2",size=1,alpha=alpha)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of pairs that have the same pathway annotation(s)",aesthetic='custom text')+
  scale_colour_manual("Distance", 
                      breaks = c("PCC based distance", "Spearman based distance", "Euclidean based distance"),
                      values = c("PCC based distance"="black", "Spearman based distance"="black", "Euclidean based distance"="black"))+
  scale_linetype_manual("Distance",values=c("PCC based distance"=2,"Spearman based distance"=3,"Euclidean based distance"=4))+
  guides(colour=guide_legend(override.aes=list(linetype=c(1,2,3)))) #ref: https://stackoverflow.com/questions/15059832/problems-with-ggplot-scale-color-and-linetype-parameters



##The corresponding plot for quantitative - "pcc" & qualitative - "MHD", "MHD collapsed condition","MI based distance"

samples=1:4000
alpha=0.5

ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,PCC_based_distance=cumsums$pcc[samples]),aes(no,PCC_based_distance,color="PCC based distance"),size=size,alpha=alpha,linetype=1)+
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative - MHD"),size=size*3,alpha=alpha,shape=16)+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative - MHD"),size=size,alpha=alpha,linetype=8)+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative collapsed condition - MHD"),size=size*3,alpha=alpha,shape=17)+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative collapsed condition - MHD"),size=size,alpha=alpha,linetype=9)+
  
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi[samples]),aes(no,CumSum,color="MI based distance"),size=size,alpha=alpha,linetype=10)+
  
  
  geom_abline(intercept=0,slope=slopeForControl,color="#D76E6E",size=1,alpha=alpha)+ #negative control
  geom_abline(intercept=0, slope=1,color="#1F89C2",size=1,alpha=alpha)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of pairs that have the same pathway annotation(s)",aesthetic='custom text')+
  scale_colour_manual("Distance", 
                      breaks = c("PCC based distance", "qualitative - MHD", "qualitative collapsed condition - MHD","MI based distance"),
                      values = c("PCC based distance"="black", "qualitative - MHD"="black", "qualitative collapsed condition - MHD"="black","MI based distance"="black"))+
  scale_linetype_manual("Distance",values=c("PCC based distance"=1,"Spearman based distance"=8,"Euclidean based distance"=9,"MI based distance"=10))+
  guides(colour=guide_legend(override.aes=list(linetype=c(1,8,9,10),shape=c(NA,16,17,NA)))) #ref: https://stackoverflow.com/questions/15059832/problems-with-ggplot-scale-color-and-linetype-parameters


#barplot at particular correlation cutoffs (for:  quantitative - "pcc" & qualitative - "MHD", "MHD collapsed condition","MI based distance"): ranking of distance=c(100,200,500,1000,2000,3000)

##Check the corresponding pcc values
cutoffs=c(100,200,500,1000,2000,3000)
sort(1-sigPheno_strain1strain2_allAnnotations_allDistances$pcc,decreasing=T)[cutoffs]


p=list()
for(i in 1:length(cutoffs)){
  dat=data.frame(Distance=c("PCC based distance", "qualitative - MHD", "qualitative collapsed condition - MHD","MI based distance","random expectation"),value=c(cumsums$pcc[cutoffs[i]],cumsums$mhd3[cutoffs[i]],cumsums$mhd3_collapsedCond[cutoffs[i]],cumsums$mi[cutoffs[i]],sum(anyCoAnnotation)/length(anyCoAnnotation)*cutoffs[i]))
  dat$Distance=factor(dat$Distance,levels=dat$Distance)#Do this to prevent automatic x axis reordering
  p[[i]]=ggplot(data=dat,aes(Distance,value))+theme_minimal()+geom_bar(stat="identity")+labs(y="Cummulative no. of co-annotations")+labs(x="")
}

gridExtra::grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],ncol=1)

dat=data.frame()
for(i in 1:length(cutoffs)){
  dat=rbind(dat,
            data.frame(Distance=c("PCC based distance", "qualitative - MHD", "qualitative collapsed condition - MHD","MI based distance","random expectation"),
                       value=c(cumsums$pcc[cutoffs[i]],cumsums$mhd3[cutoffs[i]],cumsums$mhd3_collapsedCond[cutoffs[i]],cumsums$mi[cutoffs[i]],sum(anyCoAnnotation)/length(anyCoAnnotation)*cutoffs[i]),
                       cutoff=rep(i,5)))
}
dat$Distance=factor(dat$Distance,levels=c("PCC based distance", "qualitative - MHD", "qualitative collapsed condition - MHD","MI based distance","random expectation"))#Do this to prevent automatic x axis reordering

ggplot(data=dat,aes(Distance,value,fill=cutoff))+
  geom_bar(stat="identity")+
  labs(y="Cummulative no. of co-annotations")+
  scale_fill_gradient(labels=cutoffs,low="#5F06CD",high="#E0DFF7")


