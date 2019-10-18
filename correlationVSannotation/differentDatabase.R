#Plot Corr_VS_Annot for different annotation sets

##change pheno_dat to get 2 different results (1 is using all the strains and the other is using strains that have at least 1 significant phenotype)
pheno_dat=strain1strain2_allAnnotations_allDistances
#pheno_dat=sigPheno_strain1strain2_allAnnotations_allDistances


#AnyCoAnnotation or other combinations of databases (later the experiment can be expanded to include any 2,3 databases)

CoAnnotation_list=list()
CoAnnotation_list$any=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>0,T,F)
CoAnnotation_list$Pwy=ifelse(pheno_dat[,"Pwy"]>0,T,F)
CoAnnotation_list$pcomplex=ifelse(pheno_dat[,"pcomplex"]>0,T,F)
CoAnnotation_list$operon=ifelse(pheno_dat[,"operon"]>0,T,F)
CoAnnotation_list$regulator=ifelse(pheno_dat[,"regulator"]>0,T,F)
CoAnnotation_list$kegg_modules=ifelse(pheno_dat[,"kegg_modules"]>0,T,F)
CoAnnotation_list$Pwy_and_kegg_modules=ifelse(rowSums(pheno_dat[,c("Pwy","kegg_modules")])==2,T,F)
CoAnnotation_list$Pwy_or_kegg_modules=ifelse(rowSums(pheno_dat[,c("Pwy","kegg_modules")])>0,T,F)

CoAnnotation_list$anyTwo=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>=2,T,F)
CoAnnotation_list$anyThree=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>=3,T,F)
CoAnnotation_list$anyFour=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>=4,T,F)
CoAnnotation_list$anyFive=ifelse(rowSums(pheno_dat[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])==5,T,F)




setwd("/Users/peterwu/Dropbox/Nichols_Data_mining/correlationVSannotation/figs")

#A gigantic for loop that prints every fig

start=Sys.time()
for(j in 1:length(CoAnnotation_list)){
  CoAnnotation=CoAnnotation_list[[j]]
  
  
  
  
  cumsums=list()
  distances=c("pcc","spearman","euclidean","euclidean_qualitative","euclidean_collapsedCond","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond")
  for(distance in distances){
    cumsums[[distance]]=cumsum(CoAnnotation[order(pheno_dat[,distance])])
  }
  
  str(cumsums)
  
  
  
  
  
  #plot it
  index_MHD3=pheno_dat$mhd3 %>% table %>% cumsum 
  index_collapsed_MHD3=pheno_dat$mhd3_collapsedCond %>% table %>% cumsum
  slopeForControl=sum(CoAnnotation)/length(CoAnnotation)
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  
  ##The corresponding plot for "pcc","spearman","euclidean"
  
  samples=1:4000
  alpha=0.5
  size=1
  
  
  #barplot at particular correlation cutoffs (for: "pcc","spearman","euclidean"): ranking of distance=c(100,200,500,1000,2000,3000)
  
  ##Check the corresponding pcc values
  cutoffs=c(100,200,500,1000,2000,3000)
  sort(1-pheno_dat$pcc,decreasing=T)[cutoffs]
  
  p=list()
  for(i in 1:length(cutoffs)){
    dat=data.frame(Distance=c("Pearson","Spearman","Euclidean","Random Exp."),value=c(cumsums$pcc[cutoffs[i]],cumsums$spearman[cutoffs[i]],cumsums$euclidean[cutoffs[i]],sum(CoAnnotation)/length(CoAnnotation)*cutoffs[i]))
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
  
  ggsave(paste(names(CoAnnotation_list)[j],".png",sep=""),width=12,height=8)
  
  
  
  
  
  #barplot at particular correlation cutoffs (for: "pcc","spearman","euclidean"): ranking of distance=c(100,200,500,1000,2000,3000)
  
  ##Check the corresponding pcc values
  cutoffs=c(100,200,500,1000,2000,3000)
  sort(1-pheno_dat$pcc,decreasing=T)[cutoffs]
  
  
  p=list()
  for(i in 1:length(cutoffs)){
    dat=data.frame(Distance=c("Pearson","Spearman","Euclidean","Random Exp."),value=c(cumsums$pcc[cutoffs[i]],cumsums$spearman[cutoffs[i]],cumsums$euclidean[cutoffs[i]],sum(CoAnnotation)/length(CoAnnotation)*cutoffs[i]))
    dat$Distance=factor(dat$Distance,levels=dat$Distance)#Do this to prevent automatic x axis reordering
    p[[i]]=ggplot(data=dat,aes(Distance,value,fill=Distance))+theme_minimal()+
      theme(axis.text.x = element_text(angle = 90,hjust=1,size=15),
            plot.margin = unit(c(20,0,0.5,0), "points"))+
      scale_fill_manual(values = c("#1aa548", "#bc2ec9", "#baa118","#D76E6E"))+
      geom_bar(stat="identity",show.legend=F)+
      labs(y="")+
      labs(x="")
  }
  
  p_all=gridExtra::grid.arrange(grobs=p,ncol=length(p))
  ggsave(paste(names(CoAnnotation_list)[j],"_bar.png",sep=""),plot=p_all,width=12,height=8)
  
  
  
  
  ##The corresponding plot for quantitative - "pcc" & qualitative - "MHD", "MHD collapsed condition","MI based distance"
  
  
  
  
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
  
  
  ggsave(paste(names(CoAnnotation_list)[j],"_qualitative.png",sep=""),width=12,height=8)
  
  
  
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
                           sum(CoAnnotation)/length(CoAnnotation)*cutoffs[i]))
    dat$Distance=factor(dat$Distance,levels=dat$Distance)#Do this to prevent automatic x axis reordering
    
    p[[i]]=ggplot(data=dat,aes(Distance,value,fill=Distance))+theme_minimal()+
      theme(axis.text.x = element_text(angle = 90,hjust=1,size=15),
            plot.margin = unit(c(20,0,0.5,0), "points"))+
      scale_fill_manual(values = c("#1aa548","#ce946b","#e26b16","#8b8db2","#474368","#D76E6E"))+
      geom_bar(stat="identity",show.legend=F)+
      labs(x="",y="")
  }
  
  gridExtra::grid.arrange(grobs=p,ncol=length(p))
  
  ggsave(paste(names(CoAnnotation_list)[j],"_qualitative_bar.png",sep=""),plot=p_all,width=12,height=8)
  
  
  
  
}
end=Sys.time()
end-start #Time difference of 3.470594 mins










