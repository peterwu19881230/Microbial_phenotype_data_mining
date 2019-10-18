#Prelim figs. Some of them are just copied-pasted from existing scripts

#avg |PCC| for EcoCyc pathways (originally written in sigPehno_allAnalysis.R)
##Note: Here I use just pwys in Nichols' fig.S1


#boxplot for n=c(3:22,27,28,31,41,47,48)
#Should I do this for each box?: no matter how many genes are used, 
# pathways of same No. of genes defined by EcoCyc should be put into 1 box 

#use n gene pathway according to Pathways.col -> box plot for n=1,2,....max(n) gene pathway in Pathways.col

#boxplot
#Note: up to this point there are only 280 pathways being used
ordered_cors_for_pwys=abs_cors_for_pwys[order(abs_cors_for_pwys$no_gene_used,-abs_cors_for_pwys$avg_pcc),]
##Ref about using order(a,-1): https://stackoverflow.com/questions/7793295/how-to-order-a-data-frame-by-one-descending-and-one-ascending-column
##Note: I feel that the part about rev() is BSing. rev() doesn't really help (I wonder if the person really had run it)

#have to turn no_gene_used into a factor in order to get proper color filling
ordered_cors_for_pwys$no_gene_used=as.factor(ordered_cors_for_pwys$no_gene_used)
ordered_cors_for_pwys$pwy=as.factor(ordered_cors_for_pwys$pwy)

tab=left_join(df.abs_cors_for_pwys_all,ordered_cors_for_pwys[,c("pwy","no_gene_used")],by="pwy")
tab=arrange(tab,as.numeric(no_gene_used),-abs_pcc)
tab$pwy=factor(tab$pwy,levels=unique(tab$pwy))#Do this to prevent automatic x axis reordering

#Have to figure out how to do: geom_segment(x=,xend=,y=,yend=) to each bar
ggplot(tab[tab$no_gene_used %in% c(3,4,5,6,7),], aes(x=pwy, y=abs_pcc,colour=no_gene_used,alpha=0.6)) + 
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


ggplot(tab[!tab$no_gene_used %in% c(3,4,5,6,7),], aes(x=pwy, y=abs_pcc,fill=no_gene_used)) + 
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)





#Correlation V.S. annotation with all 3979 strains (originally written in correlationVSannotation/allAnalysis.R) 
#Correlation V.S. annotation with strains that have at least 1 significant phenotypes (originally written in correlationVSannotation/sigPheno_allAnalysis.R)

##change pheno_dat to get 2 different results (1 is using all the strains and the other is using strains that have at least 1 significant phenotype)
pheno_dat=strain1strain2_allAnnotations_allDistances
#pheno_dat=sigPheno_strain1strain2_allAnnotations_allDistances



anyCoAnnotation=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>0,
                       T,F)

cumsums=list()
distances=c("pcc","spearman","euclidean","euclidean_qualitative","euclidean_collapsedCond","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond")
for(distance in distances){
  cumsums[[distance]]=cumsum(anyCoAnnotation[order(pheno_dat[,distance])])
}

str(cumsums)





#plot it
index_MHD3=pheno_dat$mhd3 %>% table %>% cumsum 
index_collapsed_MHD3=pheno_dat$mhd3_collapsedCond %>% table %>% cumsum
slopeForControl=sum(anyCoAnnotation)/length(anyCoAnnotation)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work

##The corresponding plot for "pcc","spearman","euclidean"

samples=1:4000
alpha=0.5
size=1

#grey scale
ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,PCC_based_distance=cumsums$pcc[samples]),aes(no,PCC_based_distance,color="PCC based distance"),size=size,alpha=alpha,linetype=1)+
  geom_line(data=data.frame(no=samples,Spearman_based_distance=cumsums$spearman[samples]),aes(no,Spearman_based_distance,color="Spearman based distance"),size=size,alpha=alpha,linetype=2)+
  geom_line(data=data.frame(no=samples,Euclidean_based_distance=cumsums$euclidean[samples]),aes(no,Euclidean_based_distance,color="Euclidean based distance"),size=size,alpha=alpha,linetype=3)+
  
  geom_abline(intercept=0,slope=slopeForControl,color="#D76E6E",size=1,alpha=alpha)+ #negative control
  geom_abline(intercept=0, slope=1,color="#1F89C2",size=1,alpha=alpha)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of co-annotated pairs",aesthetic='custom text')+
  scale_colour_manual("Distance", 
                      breaks = c("PCC based distance", "Spearman based distance", "Euclidean based distance"),
                      values = c("PCC based distance"="black", "Spearman based distance"="black", "Euclidean based distance"="black"))+
  scale_linetype_manual("Distance",values=c("PCC based distance"=2,"Spearman based distance"=3,"Euclidean based distance"=4))+
  guides(colour=guide_legend(override.aes=list(linetype=c(1,2,3)))) #ref: https://stackoverflow.com/questions/15059832/problems-with-ggplot-scale-color-and-linetype-parameters


#barplot at particular correlation cutoffs (for: "pcc","spearman","euclidean"): ranking of distance=c(100,200,500,1000,2000,3000)

##Check the corresponding pcc values
cutoffs=c(100,200,500,1000,2000,3000)
sort(1-pheno_dat$pcc,decreasing=T)[cutoffs]

p=list()
for(i in 1:length(cutoffs)){
  dat=data.frame(Distance=c("Pearson","Spearman","Euclidean","Random Exp."),value=c(cumsums$pcc[cutoffs[i]],cumsums$spearman[cutoffs[i]],cumsums$euclidean[cutoffs[i]],sum(anyCoAnnotation)/length(anyCoAnnotation)*cutoffs[i]))
  dat$Distance=factor(dat$Distance,levels=dat$Distance)#Do this to prevent automatic x axis reordering
  p[[i]]=ggplot(data=dat,aes(Distance,value))+theme_minimal()+
    theme(axis.text.x = element_text(angle = 90,hjust=1,size=15),
          plot.margin = unit(c(20,0,0.5,0), "points"))+
    geom_bar(stat="identity")+
    labs(y="")+
    labs(x="")
}


#colorful
ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,PCC_based_distance=cumsums$pcc[samples]),aes(no,PCC_based_distance,color="PCC based distance"),size=size,alpha=alpha)+
  geom_line(data=data.frame(no=samples,Spearman_based_distance=cumsums$spearman[samples]),aes(no,Spearman_based_distance,color="Spearman based distance"),size=size,alpha=alpha)+
  geom_line(data=data.frame(no=samples,Euclidean_based_distance=cumsums$euclidean[samples]),aes(no,Euclidean_based_distance,color="Euclidean based distance"),size=size,alpha=alpha)+
  
  geom_abline(intercept=0,slope=slopeForControl,color="#D76E6E",linetype='dashed',size=1)+ #negative control
  geom_abline(intercept=0, slope=1,color="#1F89C2",linetype='dashed',size=1)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of co-annotated pairs",aesthetic='custom text')+
  scale_colour_manual("Distance", 
                      breaks = c("PCC based distance", "Spearman based distance", "Euclidean based distance"),
                      values = c("PCC based distance"="#1aa548", "Spearman based distance"="#bc2ec9", "Euclidean based distance"="#baa118"))


#barplot at particular correlation cutoffs (for: "pcc","spearman","euclidean"): ranking of distance=c(100,200,500,1000,2000,3000)

##Check the corresponding pcc values
cutoffs=c(100,200,500,1000,2000,3000)
sort(1-pheno_dat$pcc,decreasing=T)[cutoffs]


p=list()
for(i in 1:length(cutoffs)){
  dat=data.frame(Distance=c("Pearson","Spearman","Euclidean","Random Exp."),value=c(cumsums$pcc[cutoffs[i]],cumsums$spearman[cutoffs[i]],cumsums$euclidean[cutoffs[i]],sum(anyCoAnnotation)/length(anyCoAnnotation)*cutoffs[i]))
  dat$Distance=factor(dat$Distance,levels=dat$Distance)#Do this to prevent automatic x axis reordering
  p[[i]]=ggplot(data=dat,aes(Distance,value,fill=Distance))+theme_minimal()+
    theme(axis.text.x = element_text(angle = 90,hjust=1,size=15),
          plot.margin = unit(c(20,0,0.5,0), "points"))+
    scale_fill_manual(values = c("#1aa548", "#bc2ec9", "#baa118","#D76E6E"))+
    geom_bar(stat="identity",show.legend=F)+
    labs(y="")+
    labs(x="")
}

gridExtra::grid.arrange(grobs=p,ncol=length(p))





##The corresponding plot for quantitative - "pcc" & qualitative - "MHD", "MHD collapsed condition","MI based distance"

#grey scale
samples=1:5000 #5000 is just for not letting MHD stop at the middle
alpha=0.5

ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,PCC_based_distance=cumsums$pcc[samples]),aes(no,PCC_based_distance,color="PCC based distance"),size=size,alpha=alpha,linetype=1)+
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative - MHD"),size=size*3,alpha=alpha,shape=16)+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative - MHD"),size=size,alpha=alpha,linetype=8)+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative collapsed condition - MHD"),size=size*3,alpha=alpha,shape=17)+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative collapsed condition - MHD"),size=size,alpha=alpha,linetype=9)+
  
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi[samples]),aes(no,CumSum,color="MI based distance"),size=size,alpha=alpha,linetype=10)+
  
  
  geom_abline(intercept=0,slope=slopeForControl,color="#D76E6E",size=1,alpha=alpha)+ #negative control
  geom_abline(intercept=0, slope=1,color="#1F89C2",size=1,alpha=alpha)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of co-annotated pairs",aesthetic='custom text')+
  scale_colour_manual("Distance", 
                      breaks = c("PCC based distance", "qualitative - MHD", "qualitative collapsed condition - MHD","MI based distance"),
                      values = c("PCC based distance"="black", "qualitative - MHD"="black", "qualitative collapsed condition - MHD"="black","MI based distance"="black"))+
  scale_linetype_manual("Distance",values=c("PCC based distance"=1,"Spearman based distance"=8,"Euclidean based distance"=9,"MI based distance"=10))+
  guides(colour=guide_legend(override.aes=list(linetype=c(1,8,9,10),shape=c(NA,16,17,NA)))) #ref: https://stackoverflow.com/questions/15059832/problems-with-ggplot-scale-color-and-linetype-parameters


#barplot at particular correlation cutoffs (for:  quantitative - "pcc" & qualitative - "MHD", "MHD collapsed condition","MI based distance"): ranking of distance=c(100,200,500,1000,2000,3000)

##Check the corresponding pcc values
cutoffs=c(100,200,500,1000,2000,3000)
sort(1-pheno_dat$pcc,decreasing=T)[cutoffs]


p=list()
for(i in 1:length(cutoffs)){
  dat=data.frame(Distance=c("quantitative - PCC", "MHD", "Collapsed Cond. - MHD","MI","Collapsed Cond. - MI","Random Exp."),
                 value=c(cumsums$pcc[cutoffs[i]],cumsums$mhd3[cutoffs[i]],cumsums$mhd3_collapsedCond[cutoffs[i]],
                         cumsums$mi[cutoffs[i]],cumsums$mi_collapsedCond[cutoffs[[i]]],
                         sum(anyCoAnnotation)/length(anyCoAnnotation)*cutoffs[i]))
  dat$Distance=factor(dat$Distance,levels=dat$Distance)#Do this to prevent automatic x axis reordering
  
  p[[i]]=ggplot(data=dat,aes(Distance,value))+theme_minimal()+
    theme(axis.text.x = element_text(angle = 90,hjust=1,size=15),
          plot.margin = unit(c(20,0,0.5,0), "points"))+
    geom_bar(stat="identity")+
    labs(y="")+
    labs(x="")
}

gridExtra::grid.arrange(grobs=p,ncol=length(p)) 



#colorful
samples=1:5000 #5000 is just for not letting MHD stop at the middle
alpha=0.5

ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
  geom_line(data=data.frame(no=samples,PCC_based_distance=cumsums$pcc[samples]),aes(no,PCC_based_distance,color="PCC based distance"),size=size,alpha=alpha)+
  geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative - MHD"),size=size*3,alpha=alpha,shape=16)+
  geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],CumSum=cumsums$mhd3[index_MHD3[index_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative - MHD"),size=size,alpha=alpha)+
  
  geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative collapsed condition - MHD"),size=size*3,alpha=alpha,shape=17)+
  geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],CumSum=cumsums$mhd3_collapsedCond[index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,CumSum,color="qualitative collapsed condition - MHD"),size=size,alpha=alpha)+
  
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi[samples]),aes(no,CumSum,color="MI based distance"),size=size,alpha=alpha)+
  geom_line(data=data.frame(no=samples,CumSum=cumsums$mi_collapsedCond[samples]),aes(no,CumSum,color="qualitative collapsed condition - MI"),size=size,alpha=alpha)+
  
  geom_abline(intercept=0,slope=slopeForControl,color="#D76E6E",linetype='dashed',size=1)+ #negative control
  geom_abline(intercept=0, slope=1,color="#1F89C2",linetype='dashed',size=1)+ #positive control
  scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
  labs(x = "Ranking of Distance",y="Cumulative No. of co-annotated pairs",aesthetic='custom text')+
  scale_colour_manual("Distance", 
                      breaks = c("PCC based distance", "qualitative - MHD", "qualitative collapsed condition - MHD","MI based distance","qualitative collapsed condition - MI"),
                      values = c("PCC based distance"="#1aa548","qualitative - MHD"="#ce946b",
                                 "qualitative collapsed condition - MHD"="#e26b16","MI based distance"="#8b8db2","qualitative collapsed condition - MI"="#474368"))+
  guides(colour=guide_legend(override.aes=list(shape=c(NA,16,17,NA,NA)))) #ref: https://stackoverflow.com/questions/15059832/problems-with-ggplot-scale-color-and-linetype-parameters



#barplot at particular correlation cutoffs (for: "pcc","spearman","euclidean"): ranking of distance=c(100,200,500,1000,2000,3000)

##Check the corresponding pcc values
cutoffs=c(100,200,500,1000,2000,3000)
sort(1-pheno_dat$pcc,decreasing=T)[cutoffs]


p=list()
#four_color=c("#1aa548","#bc2ec9","#baa118","#D76E6E")
for(i in 1:length(cutoffs)){
  dat=data.frame(Distance=c("quantitative - PCC", "MHD", "Collapsed Cond. - MHD","MI","Collapsed Cond. - MI","Random Exp."),
                 value=c(cumsums$pcc[cutoffs[i]],cumsums$mhd3[cutoffs[i]],cumsums$mhd3_collapsedCond[cutoffs[i]],
                         cumsums$mi[cutoffs[i]],cumsums$mi_collapsedCond[cutoffs[[i]]],
                         sum(anyCoAnnotation)/length(anyCoAnnotation)*cutoffs[i]))
  dat$Distance=factor(dat$Distance,levels=dat$Distance)#Do this to prevent automatic x axis reordering
  
  p[[i]]=ggplot(data=dat,aes(Distance,value,fill=Distance))+theme_minimal()+
    theme(axis.text.x = element_text(angle = 90,hjust=1,size=15),
          plot.margin = unit(c(20,0,0.5,0), "points"))+
    scale_fill_manual(values = c("#1aa548","#ce946b","#e26b16","#8b8db2","#474368","#D76E6E"))+
    geom_bar(stat="identity",show.legend=F)+
    labs(x="",y="")
}

gridExtra::grid.arrange(grobs=p,ncol=length(p))






