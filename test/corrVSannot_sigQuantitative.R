#Use ggplot to do correlation VS annotation

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)

source(paste(dir_of_workingScript,"/distanceForDiffTypeOfPhenotype.R",sep="")


#subset the data so only strains that have at least 1 significant phenotype are used
dat=strain1strain2_allAnnotations_allDistances[strain1strain2_allAnnotations_allDistances$sigPhenoInBoth==T,c("strain1","strain2","Pwy")]


# Use original phenotype data
pcc=strain1strain2_allAnnotations_allDistances$pcc[strain1strain2_allAnnotations_allDistances$sigPhenoInBoth==T]
cumSum_original=dat$Pwy[order(pcc)] %>% cumsum

# This is used: distances__sigQuantitative_Data_324cutoff 
Table=left_join(dat,distances__sigQuantitative_Data_324cutoff,by=c("strain1","strain2"))
cumSum_1=strain1_strain2_pcc$Pwy[order(strain1_strain2_pcc$pcc)] %>% cumsum


# This is used: distances__sigQuantitative_Data_324cutoff_condCollapsed
Table2=left_join(dat,distances__sigQuantitative_Data_324cutoff_condCollapsed,by=c("strain1","strain2"))
cumSum_2=strain1_strain2_pcc$Pwy[order(strain1_strain2_pcc$pcc)] %>% cumsum



#(!)Have to solve this problem later: ranking for the spearman, euclidean should be recoded because there are many same distances
samples=1:4000
#samples=1:length(cumSum_original)
alpha=0.5
size=1  
  
slope_negativeControl=sum(dat$Pwy)/length(dat$Pwy)

library(scales) #This is for label=comma (not using scientific notation and put commas) to work
    
plot_result=function(samples=1:4000){
  
  ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
    geom_line(data=data.frame(no=samples,cumSum_original=cumSum_original[samples]),aes(no,cumSum_original,color="Original data"),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,cumSum_1=cumSum_1[samples]),aes(no,cumSum_1,color="non-significant fitness removed"),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,cumSum_2=cumSum_2[samples]),aes(no,cumSum_2,color="non-significant fitness removed + collapse conditions"),size=size,alpha=alpha)+
    
    ##I tried to show legend for geom_abline but even the following failed
    ##https://groups.google.com/forum/#!topic/ggplot2/HtMK2ynltyw
    geom_abline(intercept=0,slope=slope_negativeControl,color="#D76E6E",linetype='dashed',size=1)+ #negative control
    geom_abline(intercept=0, slope=1,color="#1F89C2",linetype='dashed',size=1)+ #positive control
    scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
    labs(x = "Ranking of Distance",y="Cumulative No. of co-annotated pairs",aesthetic='custom text')+
    scale_colour_manual("Distance", 
                        breaks = c("cumSum_original"="Original data", "cumSum_1"="non-significant fitness removed", "cumSum_2"="non-significant fitness removed + collapse conditions"),
                        values = c("Original data"="#1aa548", "non-significant fitness removed"="#bc2ec9", "non-significant fitness removed + collapse conditions"="#baa118"))
  
  #+guides(colour=guide_legend(override.aes=list(shape=c(NA,NA,NA,NA,16,NA)))) #ref: https://stackoverflow.com/questions/15059832/problems-with-ggplot-scale-color-and-linetype-parameters
  
}

p_first_10000=plot_result(samples=1:10000)
p=plot_result(samples=1:length(cumSum_original))


#png(filename="corrVSannot_significantQuantitative_first_10000.png",height=600,width=1000)
#p_first_10000
#dev.off()

#png(filename="corrVSannot_significantQuantitative.png",height=600,width=1000) #This takes about 5 min
#p
#dev.off()

  
#explanation for the bad result:

#consider this pair in the phenotype data sigQuantitative_Data_324cutoff:  strain 5, 47, which has pcc=1
strain5=sigQuantitative_Data_324cutoff[5,] %>% as.numeric
strain47=sigQuantitative_Data_324cutoff[47,] %>% as.numeric
cor(strain5,strain47) #1


#It's pcc in All_Data_NAimputed:
strain5_original=All_Data_NAimputed[5,] %>% as.numeric
strain47_original=All_Data_NAimputed[47,] %>% as.numeric
cor(strain5_original,strain47_original) #~0.28  
    

#Create a table to compare distances in the original dataset to the pcc in truncated datasets
dat=strain1strain2_allAnnotations_allDistances[strain1strain2_allAnnotations_allDistances$sigPhenoInBoth==T,c("strain1","strain2","Pwy",
                                                                                                              "pcc","mhd3","mhd3_collapsedCond","mi","mi_collapsedCond")]

strain1_strain2_sigQuantPCC=distances__sigQuantitative_Data_324cutoff[,c("strain1","strain2","pcc")]
names(strain1_strain2_sigQuantPCC)[3]="sigQuantPCC"

strain1_strain2_sigQuantCondCollapsedPCC=distances__sigQuantitative_Data_324cutoff_condCollapsed[,c("strain1","strain2","pcc")]
names(strain1_strain2_sigQuantCondCollapsedPCC)[3]="sigQuantCondCollapsedPCC"


Table_diffPCC=left_join(dat,strain1_strain2_sigQuantPCC,by=c("strain1","strain2"))
Table_diffPCC=full_join(Table_diffPCC,strain1_strain2_sigQuantCondCollapsedPCC,by=c("strain1","strain2"))    
    

