
distance_column="mhd3"

distance=strain1strain2_allAnnotations_allDistances[[distance_column]] 


##Pwy
coAnnotated= ( strain1strain2_allAnnotations_allDistances$Pwy==1 )
pwy_distance=distance[coAnnotated] 
not_pwy_distance=distance[!coAnnotated]

##pcomplex
coAnnotated= ( strain1strain2_allAnnotations_allDistances$pcomplex==1 )
pcomplex_distance=distance[coAnnotated]
not_pcomplex_distance=distance[!coAnnotated]

##Pwy and pcomplex 
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex")])==2 )
pwyANDpcomplex_distance=distance[coAnnotated]
not_pwyANDpcomplex_distance=distance[!coAnnotated]


#all_annotSet
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5 )
all_annotSet_distance=distance[coAnnotated]
not_all_annotSet_distance=distance[!coAnnotated]


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

#bargraph for the means
barplot(c(324,322.4,322.4,319.8,318.4))



xlabs=c("All gene pairs","Same pathways","Same protein complexes","Same pathways and protein complexes",
        "Same in all 5 sets")

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
  scale_y_continuous("MHD")+
  theme(text=element_text(size=15),
        axis.text.y=element_text(size=30),
        axis.title=element_text(size=30))

  
dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
pdf(file=paste(dir_of_workingScript,file="avgMHD.pdf",sep="/"),height=8,width=15)
p
dev.off()





#Suggested by Dr. Cai, I can use Mann-Whitney U test. And I am thinking that I might want to use 1-sided test
#Ref: https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-the-differences-between-one-tailed-and-two-tailed-tests/
#(!)Should I really do the above test? my question is: is a subgroup different than the population: Ref: https://stats.stackexchange.com/questions/30562/how-to-test-whether-subgroup-mean-differs-from-overall-group-that-includes-the

## An example:
##population=c(c(10,11,12,13,14),c(1,2,3,4,5))
##wilcox.test(c(10,11,12,13,14),c(1,2,3,4,5),alternative="greater")


#These results suggest that maybe MHD is not a good metric
wilcox.test(pwy_distance,not_pwy_distance,alternative="greater") #p-value = 1
wilcox.test(pcomplex_distance,not_pcomplex_distance,alternative="greater") #p-value = 1
wilcox.test(pwyANDpcomplex_distance,not_pwyANDpcomplex_distance,alternative="greater") #p-value = 1
wilcox.test(all_annotSet_distance,not_all_annotSet_distance,alternative="greater") #p-value = 1






       

