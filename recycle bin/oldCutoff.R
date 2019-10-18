#sensitivity VS. precision (axes can be changed to have different metrics on x and y)
metricByCutoff=function(cutoffs,metric_list,pcc_dist){
  
  
  result_list=list()
  count=1
  for(i in seq(metric_list)){
    
    for(cutoff in cutoffs){
      Index=sum(pcc_dist<cutoff|pcc_dist==cutoff)
      
      result_list[[count]]=data.frame(sensitivity=metric_list[[i]]$sensitivity[Index],
                                      specificity=metric_list[[i]]$specificity[Index],
                                      precision=metric_list[[i]]$precision[Index],
                                      accuracy=metric_list[[i]]$accuracy[Index],
                                      group=names(metric_list)[i])
      
      count=count+1
    }
    
  }
  
  result=Reduce(rbind,result_list) 
  result$group=factor(result$group,levels=unique(result$group)) #this is to prevent ggplot to reorder legend
  
  return(result)
}


dist_list=list(
  "MHD"=TernaryCondCollapsed_mhd3_confusionMatrix_metrics,
  "pcc"=pcc_confusionMatrix_metrics,
  "spearman"=spearman_confusionMatrix_metrics,
  "MI"=mi_confusionMatrix_metrics,
  "Ternary MI"=TernaryCondCollapsed_mi_confusionMatrix_metrics
)



pcc_dist=sort(strain1strain2_allAnnotations_allDistances$pcc)





#Within 1 dataset the followin cutoffs can be choosen based on the same ranking of the pairs
cutoffs=1-c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1) #pcc #doing 1-c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1) explains why we choose pcc=0.6 as the cutoff


result=metricByCutoff(cutoffs=cutoffs,metric_list=dist_list,pcc_dist=pcc_dist)



Cols=c("black","#09D38A","#564c99","#832B05","#EC8D3A") 


p=ggplot(data=result,aes(sensitivity,precision,color=group))+
  theme_minimal()+
  geom_line()+
  geom_point(size=3)+
  theme(plot.title= element_text(size = 20),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=20))+
  scale_y_continuous(trans='log10',limits = c(0.00001,1))+
  scale_x_continuous(trans='log10',limits = c(0.00001,1))+
  scale_color_manual(values=Cols)+
  labs(color="Type of data - Correlation")+
  coord_flip() #this flips x and y axes



#No longer need the following
'
mhd_dist=sort(strain1strain2_allAnnotations_allDistances$mhd3)
pcc_dist=sort(strain1strain2_allAnnotations_allDistances$pcc)
spearman_dist=sort(strain1strain2_allAnnotations_allDistances$spearman)
mi_dist=sort(strain1strain2_allAnnotations_allDistances$mi)
mi_ternary_dist=sort(strain1strain2_allAnnotations_allDistances$mi_ternary)


distance_list=list(mhd_dist,pcc_dist,spearman_dist,mi_dist,mi_ternary_dist)
'

