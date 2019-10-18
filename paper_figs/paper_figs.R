#Generate figs in Nichols' reanalysis paper


#heatmap
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)

##quantitative
c(min(All_Data,na.rm=T),0,max(All_Data,na.rm=T)) #scope of the fitness scores
##=> From this I will set (max,min)=(10,-40)

tiff(paste(currentDir,"Heatmap_All_Data.tif",sep="/"),width = 480,height = 480*3979/324)

colRange=c(-40,0,10)
textRange=c(-40,-30,-20,-10,0,10)
color=c("blue","white","red")

Heatmap(
  All_Data,
  name="",
  col=colorRamp2(colRange,color),
  
  heatmap_legend_param = list(
    title = "", at = textRange
  ),
  
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=FALSE,
  show_column_names=FALSE,
  show_heatmap_legend = FALSE
)

dev.off()


tiff(paste(currentDir,"Heatmap_All_Data_first50X50.tif",sep="/")) #add res=500 to get a high-resolution bar

Heatmap(
  All_Data[1:50,1:50],
  name="",
  col=colorRamp2(colRange,color),
  
  heatmap_legend_param = list(
    title = "", at = textRange
  ),
  
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=FALSE,
  show_column_names=FALSE,
  show_heatmap_legend = FALSE #turn this to true to show the legend bar
)

dev.off()


##qualitative
tiff(paste(currentDir,"Heatmap_Ternary.tif",sep="/"),width = 480,height = 480*3979/324) 
Heatmap(
  Ternary_Data_NAnotimputed,
  name="",
  col=structure(c("red","white","blue"),names=c("1","0","-1")),
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=FALSE,
  show_column_names=FALSE,
  show_heatmap_legend = FALSE
)
dev.off()

tiff(paste(currentDir,"Heatmap_Ternary_first50X50.tif",sep="/")) #add res=500 to get a high-resolution bar
Heatmap(
  Ternary_Data_NAnotimputed[1:50,1:50],
  name="",
  col=structure(c("red","white","blue"),names=c("1","0","-1")), 
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=FALSE,
  show_column_names=FALSE,
  show_heatmap_legend = FALSE #turn this to true to show the legend bar
)
dev.off()


#boxplots




#Violin plot (also done in AvgPCC_pwy&ptcomplex.R)

#==============================================================================================================================


distance_column="pcc"



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
  scale_y_continuous("|PCC|")


pdf(file="avgPCC.pdf",height=8,width=15)
p
dev.off()
#==============================================================================================================================


#sensitivity, specificity, precision, accuracy
#This is in confusioinMatrix_analysis.R. May or May not be used.



#Correlation VS. Annotation (modified from prelim.R)
#==============================================================================================================================


#~~~~~  Some object prepareations ~~~~~

##change pheno_dat to get 2 different results (1 is using all the strains and the other is using strains that have at least 1 significant phenotype)
pheno_dat=strain1strain2_allAnnotations_allDistances

##change "CoAnnotation" to select different combinations of annotation sets
###CoAnnotation=ifelse(pheno_dat[,"Pwy"]==1,T,F)
###CoAnnotation=ifelse(pheno_dat[,"pcomplex"]==1,T,F)
CoAnnotation=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex")])==2,T,F)

###CoAnnotation=ifelse(pheno_dat[,"operon"]==1,T,F)
###CoAnnotation=ifelse(pheno_dat[,"regulator"]==1,T,F)
###CoAnnotation=ifelse(pheno_dat[,"kegg_modules"]==1,T,F)

###CoAnnotation=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5,T,F)

cumsums=list()
distances=c("pcc","spearman","euclidean","euclidean_qualitative","euclidean_collapsedCond","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond")
for(distance in distances){
  cumsums[[distance]]=cumsum(CoAnnotation[order(pheno_dat[,distance])])
}

str(cumsums)


slopeForControl=sum(CoAnnotation)/length(CoAnnotation)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work

##The corresponding plot for "pcc","spearman","euclidean"

samples=1:6000
alpha=0.5
size=1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Ranked pairs VS. cummulative sum
ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,PCC_based_distance=cumsums$pcc[samples]),aes(no,PCC_based_distance,color="PCC based distance"),size=size,alpha=alpha)+
  geom_line(data=data.frame(no=samples,Spearman_based_distance=cumsums$spearman[samples]),aes(no,Spearman_based_distance,color="Spearman based distance"),size=size,alpha=alpha)+
  geom_line(data=data.frame(no=samples,Euclidean_based_distance=cumsums$euclidean[samples]),aes(no,Euclidean_based_distance,color="Euclidean based distance"),size=size,alpha=alpha)+
  
  geom_line(data=data.frame(no=samples,euclidean_collapsedCond=cumsums$euclidean_collapsedCond[samples]),aes(no,euclidean_collapsedCond,color="qualitative collapsed condition - Euclidean"),size=size,alpha=alpha)+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative collapsed condition - MHD"),size=size*3,alpha=alpha,shape=16)+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative collapsed condition - MHD"),size=size,alpha=alpha)+
  
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi_collapsedCond[samples]),aes(no,CumSum,color="qualitative collapsed condition - MI"),size=size,alpha=alpha)+
  
  ##I tried to show legend for geom_abline but even the following failed
  ##https://groups.google.com/forum/#!topic/ggplot2/HtMK2ynltyw
  geom_abline(intercept=0,slope=slopeForControl,color="#D76E6E",linetype='dashed',size=1)+ #negative control
  geom_abline(intercept=0, slope=1,color="#1F89C2",linetype='dashed',size=1)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of co-annotated pairs",aesthetic='custom text')+
  scale_colour_manual("Distance", 
                      breaks = c("PCC based distance", "Spearman based distance", "Euclidean based distance", "qualitative collapsed condition - Euclidean","qualitative collapsed condition - MHD","qualitative collapsed condition - MI"),
                      values = c("PCC based distance"="#1aa548", "Spearman based distance"="#bc2ec9", "Euclidean based distance"="#baa118", "qualitative collapsed condition - Euclidean"="#fcf50d","qualitative collapsed condition - MHD"="#e26b16","qualitative collapsed condition - MI"="#474368"))+
  guides(colour=guide_legend(override.aes=list(shape=c(NA,NA,NA,NA,16,NA)))) #ref: https://stackoverflow.com/questions/15059832/problems-with-ggplot-scale-color-and-linetype-parameters



#==============================================================================================================================






