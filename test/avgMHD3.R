distance_column="mhd3"

distance=strain1strain2_allAnnotations_allDistances[[distance_column]] #whether "scale(center=F)" is done makes no difference for visualization

##Pwy
coAnnotated= ( strain1strain2_allAnnotations_allDistances$Pwy==1 )
pwy_distance=distance[coAnnotated] #I am converting the |PCC| based distance back to just |PCC|


##pcomplex
coAnnotated= ( strain1strain2_allAnnotations_allDistances$pcomplex==1 )
pcomplex_distance=distance[coAnnotated]

##Pwy and pcomplex 
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex")])==2 )
pwyANDpcomplex_distance=distance[coAnnotated]


#all_annotSet
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5 )
all_annotSet_distance=distance[coAnnotated]


# all
All=distance



# create dfs for ggplots
df_all=data.frame(all=All) 
df1 = data.frame(pwy=pwy_distance) 
df2 = data.frame(ptcom=pcomplex_distance)
df3=data.frame(pwyANDpcomplex=pwyANDpcomplex_distance) 
df4=data.frame(all_annotSet=all_annotSet_distance) 


summary(df_all) #Median: 324 #Mean: 324
summary(df1) #Median: 324 #Mean: 322.4
summary(df2) #Median: 324 #Mean: 322.4
summary(df3) #Median: 324 #Mean: 319.8
summary(df4) #Median: 319.5 #Mean: 318.4



xlabs=c("All gene pairs","Same pathways","Same protein complexes","Same pathways and protein complexes",
        "Intersection of 5 annotation sets")

p=ggplot() +
  geom_violin(data = df_all,aes(xlabs[1],all)) +
  geom_boxplot(data = df_all,aes(xlabs[1],all),width=0.1,outlier.shape = NA)+
  geom_violin(data = df1,aes(xlabs[2],pwy)) +
  geom_boxplot(data = df1,aes(xlabs[2],pwy),width=0.1,outlier.shape = NA)+
  geom_violin(data = df2,aes(xlabs[3],ptcom))+
  geom_boxplot(data = df2,aes(xlabs[3],ptcom),width=0.1,outlier.shape = NA)+
  geom_violin(data = df3,aes(xlabs[4],pwyANDpcomplex))+
  geom_boxplot(data = df3,aes(xlabs[4],pwyANDpcomplex),width=0.1,outlier.shape = NA)+
  geom_violin(data = df4,aes(xlabs[5],all_annotSet))+
  geom_boxplot(data = df4,aes(xlabs[5],all_annotSet),width=0.1,outlier.shape = NA)+
  scale_x_discrete("",limits=xlabs)+ 
  #I want the x axis to be empty. And if I don't use this, the order is not right
  scale_y_continuous("MHD") 

pdf(file="test/avgMHD.pdf",height=8,width=15)
p
dev.off()





