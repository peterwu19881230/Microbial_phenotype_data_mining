
distance_column="spearman"

correlation=(1-strain1strain2_allAnnotations_allDistances[[distance_column]]) #whether "scale(center=F)" is done makes no difference for visualization

#Offset=10
#correlation=( (1-strain1strain2_allAnnotations_allDistances[[distance_column]])+Offset ) %>% log(base=10) #These 2 lines also don't work

##Pwy
coAnnotated= ( strain1strain2_allAnnotations_allDistances$Pwy==1 )
pwy_correlation=correlation[coAnnotated] 
not_pwy_correlation=correlation[!coAnnotated]

##pcomplex
coAnnotated= ( strain1strain2_allAnnotations_allDistances$pcomplex==1 )
pcomplex_correlation=correlation[coAnnotated]
not_pcomplex_correlation=correlation[!coAnnotated]

##Pwy and pcomplex 
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex")])==2 )
pwyANDpcomplex_correlation=correlation[coAnnotated]
not_pwyANDpcomplex_correlation=correlation[!coAnnotated]


#all_annotSet
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5 )
all_annotSet_correlation=correlation[coAnnotated]
not_all_annotSet_correlation=correlation[!coAnnotated]


# all
All=correlation



# create dfs for ggplots
df_all=data.frame(all=All) 
df1 = data.frame(pwy=pwy_correlation) 
df2 = data.frame(ptcom=pcomplex_correlation)
df3=data.frame(pwyANDpcomplex=pwyANDpcomplex_correlation) 
df4=data.frame(all_annotSet=all_annotSet_correlation) 


summary(df_all) 
summary(df1) 
summary(df2) 
summary(df3) 
summary(df4) 



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
  scale_y_continuous("Spearman")+
  theme(text=element_text(size=15),
        axis.text.y=element_text(size=30),
        axis.title=element_text(size=30))


dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
pdf(file=paste(dir_of_workingScript,file="avgSpearman.pdf",sep="/"),height=8,width=15)
p
dev.off()




#Suggested by Dr. Cai, I can use Mann-Whitney U test. And I am thinking that I might want to use 1-sided test
#Ref: https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-the-differences-between-one-tailed-and-two-tailed-tests/
#(!)Should I really do the above test? my question is: is a subgroup different than the population: Ref: https://stats.stackexchange.com/questions/30562/how-to-test-whether-subgroup-mean-differs-from-overall-group-that-includes-the

## An example:
##population=c(c(10,11,12,13,14),c(1,2,3,4,5))
##wilcox.test(c(10,11,12,13,14),c(1,2,3,4,5),alternative="greater")

wilcox.test(pwy_correlation,not_pwy_correlation,alternative="greater") #p-value < 2.2e-16
wilcox.test(pcomplex_correlation,not_pcomplex_correlation,alternative="greater") #p-value < 2.2e-16
wilcox.test(pwyANDpcomplex_correlation,not_pwyANDpcomplex_correlation,alternative="greater") #p-value < 2.2e-16
wilcox.test(all_annotSet_correlation,not_all_annotSet_correlation,alternative="greater") #p-value < 2.2e-16







