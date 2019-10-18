distance_column="mhd3"


##Pwy
coAnnotated= ( strain1strain2_allAnnotations_allDistances$Pwy==1 )
pwy_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated] #I am converting the |PCC| based distance back to just |PCC|


##pcomplex
coAnnotated= ( strain1strain2_allAnnotations_allDistances$pcomplex==1 )
pcomplex_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]

##Pwy and pcomplex 
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex")])==2 )
pwyANDpcomplex_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]


#all_annotSet
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5 )
all_annotSet_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]


# all
All=1-strain1strain2_allAnnotations_allDistances[[distance_column]]



# create dfs for ggplots
df1 = data.frame(pwy=pwy_abs_pcc) %>% scale #No difference whether I scale or not
df2 = data.frame(ptcom=pcomplex_abs_pcc) %>% scale
df3=data.frame(pwyANDpcomplex=pwyANDpcomplex_abs_pcc) %>% scale
df4=data.frame(all_annotSet=all_annotSet_abs_pcc) %>% scale
df_all=data.frame(all=All) %>% scale

summary(df1)
summary(df2)
summary(df3)
summary(df4)


xlabs=c("Gene pairs in the same pathways","Gene pairs in the same protein complexes","Gene pairs in the same pathway and pcomplex",
        "Gene pairs co-annotated in all 5 annotation sets","All gene pairs")

ggplot() +
  geom_violin(data = df1,aes(xlabs[1],pwy)) +
  geom_boxplot(data = df1,aes(xlabs[1],pwy),width=0.1,outlier.shape = NA)+
  geom_violin(data = df2,aes(xlabs[2],ptcom))+
  geom_boxplot(data = df2,aes(xlabs[2],ptcom),width=0.1,outlier.shape = NA)+
  geom_violin(data = df3,aes(xlabs[3],pwyANDpcomplex))+
  geom_boxplot(data = df3,aes(xlabs[3],pwyANDpcomplex),width=0.1,outlier.shape = NA)+
  geom_violin(data = df4,aes(xlabs[4],all_annotSet))+
  geom_boxplot(data = df4,aes(xlabs[4],all_annotSet),width=0.1,outlier.shape = NA)+
  geom_violin(data = df_all,aes(xlabs[5],all)) +
  geom_boxplot(data = df_all,aes(xlabs[5],all),width=0.1,outlier.shape = NA)+
  scale_x_discrete("",limits=xlabs)+ 
  #I want the x axis to be empty. And if I don't use this, the order is not right
  scale_y_continuous("MHD")
