#Construct avg c("pcc","spearman","euclidean","euclidean_qualitative","euclidean_collapsedCond","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond") for: 
# c("Pwy","pcomplex","operon","regulator","kegg_modules")

names(strain1strain2_allAnnotations_allDistances) 

annotations=c("Pwy","pcomplex","operon","regulator","kegg_modules")
distances=c("pcc","spearman","euclidean","euclidean_qualitative","euclidean_collapsedCond","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond")


#Avg distances in any coannotations
anyCoAnnotation=ifelse(rowSums(strain1strain2_allAnnotations_allDistances[,annotations])>0,
                       T,F) 
length(anyCoAnnotation)
sum(anyCoAnnotation)
sum(anyCoAnnotation)/length(anyCoAnnotation)


distance_coannotation=list()
for(distance in distances){
  distance_coannotation[[distance]]=strain1strain2_allAnnotations_allDistances[,distance][anyCoAnnotation]
  distance_coannotation[[paste("all_",distance,sep="")]]=strain1strain2_allAnnotations_allDistances[,distance]
}

sapply(distance_coannotation,mean) #Looks like coannotated pairs don't have significantly lower distances (higher correlations)

p=list()

  p$pcc=ggplot() +
    geom_violin(data = data.frame(1-distance_coannotation$pcc),aes("Coannotated pairs",1-distance_coannotation$pcc)) +
    geom_boxplot(data = data.frame(1-distance_coannotation$pcc),aes("Coannotated pairs",1-distance_coannotation$pcc),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(1-distance_coannotation$all_pcc),aes("All",1-distance_coannotation$all_pcc)) +
    geom_boxplot(data = data.frame(1-distance_coannotation$all_pcc),aes("All",1-distance_coannotation$all_pcc),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("PCC")

  p$spearman=ggplot() +
    geom_violin(data = data.frame(1-distance_coannotation$spearman),aes("Coannotated pairs",1-distance_coannotation$spearman)) +
    geom_boxplot(data = data.frame(1-distance_coannotation$spearman),aes("Coannotated pairs",1-distance_coannotation$spearman),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(1-distance_coannotation$all_spearman),aes("All",1-distance_coannotation$all_spearman)) +
    geom_boxplot(data = data.frame(1-distance_coannotation$all_spearman),aes("All",1-distance_coannotation$all_spearman),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("SRCC")
  
  
  p$euclidean=ggplot() +
    geom_violin(data = data.frame(distance_coannotation$euclidean),aes("Coannotated pairs",distance_coannotation$euclidean)) +
    geom_boxplot(data = data.frame(distance_coannotation$euclidean),aes("Coannotated pairs",distance_coannotation$euclidean),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(distance_coannotation$all_euclidean),aes("All",distance_coannotation$all_euclidean)) +
    geom_boxplot(data = data.frame(distance_coannotation$all_euclidean),aes("All",distance_coannotation$all_euclidean),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("Euclidean")
  
  p$euclidean_qualitative=ggplot() +
    geom_violin(data = data.frame(distance_coannotation$euclidean_qualitative),aes("Coannotated pairs",distance_coannotation$euclidean_qualitative)) +
    geom_boxplot(data = data.frame(distance_coannotation$euclidean_qualitative),aes("Coannotated pairs",distance_coannotation$euclidean_qualitative),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(distance_coannotation$all_euclidean_qualitative),aes("All",distance_coannotation$all_euclidean_qualitative)) +
    geom_boxplot(data = data.frame(distance_coannotation$all_euclidean_qualitative),aes("All",distance_coannotation$all_euclidean_qualitative),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("Euclidean Qualitative")
  
  p$euclidean_collapsedCond=ggplot() +
    geom_violin(data = data.frame(distance_coannotation$euclidean_collapsedCond),aes("Coannotated pairs",distance_coannotation$euclidean_collapsedCond)) +
    geom_boxplot(data = data.frame(distance_coannotation$euclidean_collapsedCond),aes("Coannotated pairs",distance_coannotation$euclidean_collapsedCond),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(distance_coannotation$all_euclidean_collapsedCond),aes("All",distance_coannotation$all_euclidean_collapsedCond)) +
    geom_boxplot(data = data.frame(distance_coannotation$all_euclidean_collapsedCond),aes("All",distance_coannotation$all_euclidean_collapsedCond),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("euclidean collapsedCond")
  
  
  p$mhd3=ggplot() +
    geom_violin(data = data.frame(324-distance_coannotation$mhd3),aes("Coannotated pairs",324-distance_coannotation$mhd3)) +
    geom_boxplot(data = data.frame(324-distance_coannotation$mhd3),aes("Coannotated pairs",324-distance_coannotation$mhd3),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(324-distance_coannotation$all_mhd3),aes("All",324-distance_coannotation$all_mhd3)) +
    geom_boxplot(data = data.frame(324-distance_coannotation$all_mhd3),aes("All",324-distance_coannotation$all_mhd3),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("MHD3") #If I can make this into a kind of correlation, it would be easier to compare with some other correlations
  
  p$mhd3_collapsedCond=ggplot() +
    geom_violin(data = data.frame(324-distance_coannotation$mhd3_collapsedCond),aes("Coannotated pairs",324-distance_coannotation$mhd3_collapsedCond)) +
    geom_boxplot(data = data.frame(324-distance_coannotation$mhd3_collapsedCond),aes("Coannotated pairs",324-distance_coannotation$mhd3_collapsedCond),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(324-distance_coannotation$all_mhd3_collapsedCond),aes("All",324-distance_coannotation$all_mhd3_collapsedCond)) +
    geom_boxplot(data = data.frame(324-distance_coannotation$all_mhd3_collapsedCond),aes("All",324-distance_coannotation$all_mhd3_collapsedCond),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("MHD3") #If I can make this into a kind of correlation, it would be easier to compare with some other correlations
  
  p$mi=ggplot() +
    geom_violin(data = data.frame(1-distance_coannotation$mi),aes("Coannotated pairs",1-distance_coannotation$mi)) +
    geom_boxplot(data = data.frame(1-distance_coannotation$mi),aes("Coannotated pairs",1-distance_coannotation$mi),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(1-distance_coannotation$all_mi),aes("All",1-distance_coannotation$all_mi)) +
    geom_boxplot(data = data.frame(1-distance_coannotation$all_mi),aes("All",1-distance_coannotation$all_mi),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("MI")
  
  p$mi_collapsedCond=ggplot() +
    geom_violin(data = data.frame(1-distance_coannotation$mi_collapsedCond),aes("Coannotated pairs",1-distance_coannotation$mi_collapsedCond)) +
    geom_boxplot(data = data.frame(1-distance_coannotation$mi_collapsedCond),aes("Coannotated pairs",1-distance_coannotation$mi_collapsedCond),width=0.1,outlier.shape = NA)+
    geom_violin(data = data.frame(1-distance_coannotation$all_mi_collapsedCond),aes("All",1-distance_coannotation$all_mi_collapsedCond)) +
    geom_boxplot(data = data.frame(1-distance_coannotation$all_mi_collapsedCond),aes("All",1-distance_coannotation$all_mi_collapsedCond),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=c("Coannotated pairs","All"))+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("MI Collapsed Condition")
  
  #The following code runs very slow
  #gridExtra::grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],ncol=3)



#Avg distances in every coannotations
#for(distance in distances){
#  for(annotation in annotations){
    
#  }
  
#  strain1strain2_allAnnotations_allDistances[,annotation]
#}

  
#Conclusion: By using all coannotations, there doesn't seem to be significant increase in correlations (all distance metrics show similar result)
  
  
  