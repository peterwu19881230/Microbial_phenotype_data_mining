#Goal: Determine the significance of avg pcc of pwys and ptcomplexes. 
#(!) I haven't dealt with the duplicated strain problem 
#(eg. 2 strains that have a identical name might be viewed as 2 different genes in the same pwy or protein complex)
#=>But it's just a small amount => I guess I can neglect those (Otherwise my code is going to be ugly and non-systematic)


##Violin plot (The shape is a little different than using vioplot() but overall they look similar. Is it because ggplot() gives higher resolution?)
#=========================================================================================================
distance_column="pcc"

##Pwy
coAnnotated= ( strain1strain2_allAnnotations_allDistances$Pwy==1 )
pwy_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated] #I am converting the |PCC| based distance back to just |PCC|
not_pwy_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][!coAnnotated]

##pcomplex
coAnnotated= ( strain1strain2_allAnnotations_allDistances$pcomplex==1 )
pcomplex_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]
not_pcomplex_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][!coAnnotated]


##Pwy and pcomplex 
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex")])==2 )
pwyANDpcomplex_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]
not_pwyANDpcomplex_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][!coAnnotated]


#all_annotSet
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5 )
all_annotSet_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][coAnnotated]
not_all_annotSet_abs_pcc=1-strain1strain2_allAnnotations_allDistances[[distance_column]][!coAnnotated]


# all
All=1-strain1strain2_allAnnotations_allDistances[[distance_column]]



# create dfs for ggplots
df_all=data.frame(all=All)
df1 = data.frame(pwy=pwy_abs_pcc)
df2 = data.frame(ptcom=pcomplex_abs_pcc)
df3=data.frame(pwyANDpcomplex=pwyANDpcomplex_abs_pcc)
df4=data.frame(all_annotSet=all_annotSet_abs_pcc)

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
  scale_y_continuous("|PCC|")+
  theme(text=element_text(size=15), #if this is 25 some labels will overlap each other
        axis.text.y=element_text(size=25),
        axis.title=element_text(size=25))


dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
pdf(file=paste(dir_of_workingScript,"avgPCC.pdf",sep="/"),height=8,width=15)
p
dev.off()
#An alternative to save the plot
#ggsave(paste(dir_of_workingScript,"avgPCC.pdf",sep="/"),height=8,width=15)



#Suggested by Dr. Cai, I can use Mann-Whitney U test. And I am thinking that I might want to use 1-sided test
#Ref: https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-the-differences-between-one-tailed-and-two-tailed-tests/
#(!)Should I really do the above test? my question is: is a subgroup different than the population: Ref: https://stats.stackexchange.com/questions/30562/how-to-test-whether-subgroup-mean-differs-from-overall-group-that-includes-the

## An example:
##population=c(c(10,11,12,13,14),c(1,2,3,4,5))
##wilcox.test(c(10,11,12,13,14),c(1,2,3,4,5),alternative="greater")

wilcox.test(pwy_abs_pcc,not_pwy_abs_pcc,alternative="greater") #p-value < 2.2e-16
wilcox.test(pcomplex_abs_pcc,not_pcomplex_abs_pcc,alternative="greater") #p-value < 2.2e-16
wilcox.test(pwyANDpcomplex_abs_pcc,not_pwyANDpcomplex_abs_pcc,alternative="greater") #p-value < 2.2e-16
wilcox.test(all_annotSet_abs_pcc,not_all_annotSet_abs_pcc,alternative="greater") #p-value < 2.2e-16

hist(All,breaks=100,freq=F)
hist(pwy_abs_pcc,breaks=100,freq=F,col=rgb(0,0,1,0.2),add=T)


#an alternative for violin plot
plot(density(All))
lines(density(pwy_abs_pcc))
lines(density(pcomplex_abs_pcc))
lines(density(pwyANDpcomplex_abs_pcc))
lines(density(all_annotSet_abs_pcc))

#=========================================================================================================


#load many objects from this .RData (defined in AvgPCC_pwy&ptcomplex_obj_prep.R):
## These objects include: cors_for_pwys,cors_for_pwys_all,abs_cors_for_pwys,abs_cors_for_pwys_all,cors_for_pcomplexes,cors_for_pcomplexes_all,abs_cors_for_pcomplexes,abs_cors_for_pcomplexes_all,samePwy_df,samePcomplex_df,samePwyORPcomplex_df,samePwyPcomplex_df
load("Data/cors_for_Pwy_Pcomplex.RData")


#For distribution
#======================================================================
#Distribution of mean Pwy |pcc|
hist(abs(samePwy_df$pcc),main="Distribution of intra-pathway mean |pcc|",xlab="Pearson Correlation Coefficient")
summary(abs(samePwy_df$pcc))

#Distribution of mean Pcomplex |pcc|
hist(abs(samePcomplex_df$pcc),main="Distribution of intra-protein complex mean pcc",xlab="Pearson Correlation Coefficient")
summary(abs(samePcomplex_df$pcc))



#Distribution of mean Pwy & Pcomplex |pcc| (pairs that are co-annotated by both annotation sets)
hist(abs(samePwyPcomplex_df$pcc),main="Distribution of intra-protein complex mean pcc",xlab="Pearson Correlation Coefficient")
summary(abs(samePwyPcomplex_df$pcc))

##Verify with the newly created strain1strain2_allAnnotations_allDistances 
coAnnotated= ( rowSums(strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex")])==2 )
abs_pcc=strain1strain2_allAnnotations_allDistances$pcc[coAnnotated]

length(1-abs_pcc)
summary(1-abs_pcc)

length(abs(samePwyPcomplex_df$pcc))
summary(abs(samePwyPcomplex_df$pcc))

hist(1-abs_pcc)
hist(abs(samePwyPcomplex_df$pcc))
##The histogram look similar but not exactly the same. I suspect it's because I used the original data VS. data with NA imputed







#Negative control (average of all |pcc|)
mean(abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient)) ##0.09278531



##Info of : (1) no. of pairs (2) median


lapply(list(df1,df2,df3,df4,df5),FUN = function(df){
  no=dim(df)[1]
  med=median(df[[1]]) #subset this dataframe with[[1]] gets a numeric vector (a dataframe is a special case of list, where columns are elements in the list)
  return(list(no=no,med=med))
  })


mean(abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient))




#Distribution of mean pccs
hist(samePwy_df$pcc,main="Distribution of mean pccs of 290 pathways",xlab="Pearson Correlation Coefficient")
summary(samePwy_df$pcc)


#======================================================================


#ggplot boxplot

#Random line
random_expectation=mean(1-strain1strain2_allAnnotations_allDistances$pcc)

#Function to generate df for the boxplot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pre_process_for_ggplot=function(abs_cors_for_annotations,abs_cors_for_annotations_all){
  #have to turn no_gene_used into a factor in order to get proper color filling
  abs_cors_for_annotations$no_gene_used=as.factor(abs_cors_for_annotations$no_gene_used)
  abs_cors_for_annotations[[1]]=as.factor(abs_cors_for_annotations[[1]]) #suppose the annotation column is the 1st column
  
  #Lost the way I plot using different colors for genes in n gene pathways, but I think it's not important anymore
  df.abs_cors_for_annotations_all=data.frame(annotation=rep(names(abs_cors_for_annotations_all),sapply(abs_cors_for_annotations_all,length)),abs_pcc=unlist(abs_cors_for_annotations_all))
  names(df.abs_cors_for_annotations_all)[1]=names(abs_cors_for_annotations)[1]
  
  tab=left_join(df.abs_cors_for_annotations_all,abs_cors_for_annotations[,c(names(abs_cors_for_annotations)[1],"median_pcc","no_gene_used")],by=names(abs_cors_for_annotations)[1])
  tab=arrange(tab,as.numeric(no_gene_used),-median_pcc,-abs_pcc) #tab[order(a,-b,-c),] should also work
  ##Ref about using order(a,-1): https://stackoverflow.com/questions/7793295/how-to-order-a-data-frame-by-one-descending-and-one-ascending-column
  ##Note: I feel that the part about rev() is BSing. rev() doesn't really help (I wonder if the person really had run it)
  
  
  tab[[names(abs_cors_for_annotations)[1]]]=factor(tab[[names(abs_cors_for_annotations)[1]]],levels=unique(tab[[names(abs_cors_for_annotations)[1]]]))#Do this to prevent automatic x axis reordering
  
  ##ref: https://stackoverflow.com/questions/43877663/order-multiple-variables-in-ggplot2
  ##ref: http://www.sthda.com/english/wiki/ggplot2-facet-split-a-plot-into-a-matrix-of-panels
  ##Also, I am separating them into 2 lines
  
  return(tab)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Pwy

#no. of genes in pathways= c(3:22,27,28,31,41,47,48)
#I am currently using this: no matter how many genes are used, 
# pathways of same No. of genes defined by EcoCyc should be put into 1 box 
#=========================================================================================================

#Note: up to this point there are 282 pathways being used

pwy_tab=pre_process_for_ggplot(abs_cors_for_pwys,abs_cors_for_pwys_all)
tab=pwy_tab


##for n=2,3,4,5
p_Pwy1=ggplot(tab[tab$no_gene_used %in% c(2,3,4,5),], aes(x=Pwy, y=abs_pcc)) +  
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add fill=no_gene_used to aes() to make it colorful


#I need this to correct n=2 into a dot plot
p_Pwy1_forCorrection_1=ggplot(tab[tab$no_gene_used %in% c(2,3,4,5),], aes(x=Pwy, y=abs_pcc))+ #I don't understand why "fill=no_gene_used" will give all black dots
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point(color="grey")+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful


#I need this to correct n=3 into a dot plot
p_Pwy1_forCorrection_2=ggplot(tab[tab$no_gene_used %in% c(2,3,4,5),], aes(x=Pwy, y=abs_pcc))+ #I don't understand why "fill=no_gene_used" will give all black dots
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point(color="grey")+ylab("")+xlab("")+
  stat_summary(fun.y = median,fun.ymin=median,fun.ymax = median,color="black",
               geom = "crossbar", size = 0.3,width=1)+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful



##for n=6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 26, 27, 28, 31, 41, 48
p_Pwy2=ggplot(tab[tab$no_gene_used %in% c(6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,26,27,28,31,41,48),], 
       aes(x=Pwy, y=abs_pcc)) +
  theme_minimal()+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add fill=no_gene_used to aes() to make it colorful


#=========================================================================================================



#Pcomplex
#=========================================================================================================

pcomplex_tab=pre_process_for_ggplot(abs_cors_for_pcomplexes,abs_cors_for_pcomplexes_all)
tab=pcomplex_tab

tab$no_gene_used %>% unique # 2  3  4  5  6  7  9  10 11 12 27 


#I need this to correct n=2 into a dot plot
p_Pcomplex_forCorrection_1=ggplot(tab[tab$no_gene_used %in% c(2,3),], aes(x=pcomplex, y=abs_pcc))+ #I don't understand why "fill=no_gene_used" will give all black dots
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Protein Complex")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point(color="grey")+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Protein Complex Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful


#I need this to correct n=3 into a dot plot
p_Pcomplex_forCorrection_2=ggplot(tab[tab$no_gene_used %in% c(2,3),], aes(x=pcomplex, y=abs_pcc))+ #I don't understand why "fill=no_gene_used" will give all black dots
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Protein Complex")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point(color="grey")+ylab("")+xlab("")+
  stat_summary(fun.y = median,fun.ymin=median,fun.ymax = median,color="black",
               geom = "crossbar", size = 0.3,width=1)+
  #ylab("|PCC|")+xlab("EcoCyc Protein Complex Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful



p_Pcomplex=ggplot(tab[tab$no_gene_used %in% c(4, 5, 6,  7,  9, 10, 11, 12, 27),], aes(x=pcomplex, y=abs_pcc)) + 
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Protein Complex")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Protein Complex Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add fill=no_gene_used to aes() to make it colorful


#=========================================================================================================

#Save the graph in pdf

#Get the directory of current working script (Used when in Rstudio). Ref: https://stackoverflow.com/questions/47044068/get-the-path-of-current-script/47045368
'
dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path) 

pdf(file=paste(dir_of_workingScript,"/pwy_pcom_avgPCC.pdf",sep=""),width=15,height=5) 
p_Pwy1
p_Pwy1_forCorrection_1
p_Pwy1_forCorrection_2
p_Pwy2
p_Pcomplex_forCorrection_1
p_Pcomplex_forCorrection_2
p_Pcomplex
dev.off()
'





