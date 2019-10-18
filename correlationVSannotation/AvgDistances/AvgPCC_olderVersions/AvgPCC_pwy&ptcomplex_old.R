#Goal: Determine the significance of avg pcc of pwys and ptcomplexes. 
#(!) I haven't dealt with the duplicated strain problem 
#(eg. 2 strains that have a identical name might be viewed as 2 different genes in the same pwy or protein complex)
#=>But it's just a small amount => I guess I can neglect those (Otherwise my code is going to be ugly and non-systematic)

#(Note)Same graph using different code is in paper_figs.R


#load many objects from this .RData (defined in AvgPCC_pwy&ptcomplex_obj_prep.R):
## These objects include: cors_for_pwys,cors_for_pwys_all,abs_cors_for_pwys,abs_cors_for_pwys_all,cors_for_pcomplexes,cors_for_pcomplexes_all,abs_cors_for_pcomplexes,abs_cors_for_pcomplexes_all,samePwy_df,samePcomplex_df,samePwyORPcomplex_df,samePwyPcomplex_df
load("Data/cors_for_Pwy_Pcomplex.RData")


##Violin plot (The shape is a little different than using vioplot() but overall they look similar. Is it because ggplot() gives higher resolution?)
#=========================================================================================================
### pwy
df1 = data.frame(pwy=abs(samePwy_df$pcc))

### pcomplex
df2 = data.frame(ptcom=abs(samePcomplex_df$pcc))

###  pwy "or" ptcomplex
df3=data.frame(pwyORpcomplex=abs(samePwyORPcomplex_df$pcc))

###  pwy "and" ptcomplex
df4=data.frame(pwyANDpcomplex=abs(samePwyPcomplex_df$pcc))

### all
df5=data.frame(all=abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient))


ggplot() +
  geom_violin(data = df5,aes("All gene pairs",all)) +
  geom_boxplot(data = df5,aes("All gene pairs",all),width=0.1,outlier.shape = NA)+
  geom_violin(data = df1,aes("Same pathways",pwy)) + #This means: Gene pairs in the same pathways
  geom_boxplot(data = df1,aes("Same pathways",pwy),width=0.1,outlier.shape = NA)+
  geom_violin(data = df2,aes("Same protein complexes",ptcom))+
  geom_boxplot(data = df2,aes("Same protein complexes",ptcom),width=0.1,outlier.shape = NA)+
  geom_violin(data = df4,aes("Same pathway and protein complexes",pwyANDpcomplex))+
  geom_boxplot(data = df4,aes("Same pathway and protein complexes",pwyANDpcomplex),width=0.1,outlier.shape = NA)+
  scale_x_discrete("",limits=c("All gene pairs","Same pathways","Same protein complexes","Same pathway and protein complexes"))+ 
  #I want the x axis to be empty. And if I don't use this, the order is not right
  scale_y_continuous("|PCC|")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 15),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))
#=========================================================================================================


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
hist(cor_in_pwy,main="Distribution of mean pccs of 290 pathways",xlab="Pearson Correlation Coefficient")
summary(cor_in_pwy)


#Suggested by Dr. Cai, I can use Mann-Whitney U test. And I am thinking that I might want to use 1-sided test
#Ref: https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-the-differences-between-one-tailed-and-two-tailed-tests/
#(!)Should I really do the above test? my question is: is a subgroup different than the population: https://stats.stackexchange.com/questions/30562/how-to-test-whether-subgroup-mean-differs-from-overall-group-that-includes-the

## An example:
##wilcox.test(c(10,10,10,10,10),c(1,2,3,4,5),alternative="greater")

wilcox.test(abs(cor_in_pwy),abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient),alternative="greater") #p-value < 2.2e-16
wilcox.test(abs(cor_in_pcomplex),abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient),alternative="greater") #p-value < 2.2e-16

#======================================================================


#ggplot boxplot

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


##for n=3,4,5,6,7
p_Pwy1=ggplot(tab[tab$no_gene_used %in% c(3,4,5,6,7),], aes(x=Pwy, y=abs_pcc,fill=no_gene_used)) + 
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)


#I need this to correct n=3 into a dot plot
p_Pwy1_forCorrection=ggplot(tab[tab$no_gene_used %in% c(3,4,5,6,7),], aes(x=Pwy, y=abs_pcc,colour=no_gene_used,alpha=0.6))+ #I don't understand why "fill=no_gene_used" will give all black dots
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point()+ylab("")+xlab("")+
  stat_summary(fun.y = median,fun.ymin=median,fun.ymax = median,color="black",
               geom = "crossbar", size = 0.3,width=1)+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)


##for n=8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 26, 27, 28, 31, 41, 48
p_Pwy2=ggplot(tab[tab$no_gene_used %in% c(8,9,10,11,12,13,14,15,16,17,18,20,21,22,26,27,28,31,41,48),], 
       aes(x=Pwy, y=abs_pcc,fill=no_gene_used)) + 
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)



#=========================================================================================================



#Pcomplex
#=========================================================================================================

pcomplex_tab=pre_process_for_ggplot(abs_cors_for_pcomplexes,abs_cors_for_pcomplexes_all)
tab=pcomplex_tab

tab$no_gene_used %>% unique # 3  4  5  6  7  9  10 11 12 27 #=> seems like I don't have to subset into upper and lower panels

p_Pcomplex=ggplot(tab, aes(x=pcomplex, y=abs_pcc,fill=no_gene_used)) + 
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Protein Complex")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Protein Complex Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)

#I need this to correct n=3 into a dot plot
p_Pcomplex_forCorrection=ggplot(tab, aes(x=pcomplex, y=abs_pcc,colour=no_gene_used,alpha=0.6))+ #I don't understand why "fill=no_gene_used" will give all black dots
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Protein Complex")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point()+ylab("")+xlab("")+
  stat_summary(fun.y = median,fun.ymin=median,fun.ymax = median,color="black",
               geom = "crossbar", size = 0.3,width=1)+
  #ylab("|PCC|")+xlab("EcoCyc Protein Complex Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)

#=========================================================================================================

#Save the graph in pdf

#Get the directory of current working script (Used when in Rstudio). Ref: https://stackoverflow.com/questions/47044068/get-the-path-of-current-script/47045368
'
dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path) 

pdf(file=paste(dir_of_workingScript,"/pwy_pcom_avgPCC.pdf",sep=""),width=15,height=5) 
p_Pwy1
p_Pwy1_forCorrection
p_Pwy2
p_Pcomplex
p_Pcomplex_forCorrection
dev.off()
'





