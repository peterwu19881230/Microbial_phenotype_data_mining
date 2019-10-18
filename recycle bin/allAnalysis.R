#Correlation VS. annotation from all annotation sources and by all distances

head(strain1strain2_allAnnotations_allDistances)

anyCoAnnotation=ifelse(rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>0,
                      T,F)

cumsums=list()
distances=c("pcc","spearman","euclidean","euclidean_qualitative","euclidean_collapsedCond","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond")
for(distance in distances){
  cumsums[[distance]]=cumsum(anyCoAnnotation[order(strain1strain2_allAnnotations_allDistances[,distance])])
}

str(cumsums)


#plot it
index_MHD3=strain1strain2_allAnnotations_allDistances$mhd3 %>% table %>% cumsum 
index_collapsed_MHD3=strain1strain2_allAnnotations_allDistances$mhd3_collapsedCond %>% table %>% cumsum
slopeForControl=sum(anyCoAnnotation)/length(anyCoAnnotation)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work


#samples=seq(from=0,to=7914231,by=1989); samples[1]=1
samples=1:4000
ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,CumSum=cumsums$pcc[samples]),aes(no,CumSum,color="PCC based distance"))+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$spearman[samples]),aes(no,CumSum,color="Spearman based distance"))+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$euclidean[samples]),aes(no,CumSum,color="Euclidean based distance"))+
  
  geom_line(data=data.frame(no=samples,CumSum=cumsums$euclidean_qualitative[samples]),aes(no,CumSum,color="Euclidean qualitative data"))+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$euclidean_collapsedCond[samples]),aes(no,CumSum,color="Euclidean collapsed condition"))+
  
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"))+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD"))+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"))+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="MHD collapsed condition"))+
  
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi[samples]),aes(no,CumSum,color="MI based distance"))+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi_collapsedCond[samples]),aes(no,CumSum,color="MI based distance - collapsed condition"))+
  
  geom_abline(intercept=0,slope=slopeForControl,linetype='dotted',color="#D76E6E",size=1)+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="#1F89C2",size=1)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of pairs that have the same pathway annotation(s)",aesthetic='custom text')+
  scale_color_discrete(name="Distance")





