#x: ranked pairs , y: Sensitivity --- different combinatio of annotation sets on the same graph
#(From confusionMatrix_analysis.R I noticed that sensitivuty is probably the only thing that is worth plotting)



#4 metrics defined in Wikipedia:
##sensitivity = TP / (TP +FN)
##specificity = TN / (TN+FP)
##precision = TP/(TP+FP)
##accuracy = (TP + TN)/total


##suppose the cutoff pcc=rho

##TP=no. of pairs that: (a) have pcc >= rho (b) are co-annotated
##TN=no. of pairs that: (a) have pcc < rho (b) are not co-annotated
##FP=no. of pairs that: (a) have pcc >= rho (b) are not co-annotated
##FN=no. of pairs that: (a) have pcc < rho (b) are co-annotated








#Function to get the graph.
#========================================================================================================================
plot_metrics=function(samples,alpha,size,three_dist_list,metric,
                      index_MHD3=NULL,index_collapsed_MHD3=NULL, #These are to give indices where rankings make real sense
                      x_lab="Low distance -- ranked pairs -- High distance",y,legend_title="Types of data - distance"){
  
  random_metric=paste("random",metric,sep="_")
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[1]][[metric]][samples]),aes(no,Metric,color=names(three_dist_list[1])),size=size,alpha=alpha)+
    
    geom_point(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],Metric=three_dist_list[[2]][[metric]][index_MHD3[index_MHD3<=max(samples)]]),aes(no,Metric,color=names(three_dist_list[2])),size=size*3,alpha=alpha,shape=16)+
    geom_line(data=data.frame(no=index_MHD3[index_MHD3<=max(samples)],Metric=three_dist_list[[2]][[metric]][index_MHD3[index_MHD3<=max(samples)]]),aes(no,Metric,color=names(three_dist_list[2])),size=size,alpha=alpha)+
    geom_point(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],Metric=three_dist_list[[3]][[metric]][index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,Metric,color=names(three_dist_list[3])),size=size*3,alpha=alpha,shape=16)+
    geom_line(data=data.frame(no=index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)],Metric=three_dist_list[[3]][[metric]][index_collapsed_MHD3[index_collapsed_MHD3<=max(samples)]]),aes(no,Metric,color=names(three_dist_list[3])),size=size,alpha=alpha)+
    
    
    #This is random control:
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[1]][[random_metric]][samples]),aes(no,Metric,color=names(three_dist_list[1])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[2]][[random_metric]][samples]),aes(no,Metric,color=names(three_dist_list[2])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    geom_line(data=data.frame(no=samples,Metric=three_dist_list[[3]][[random_metric]][samples]),aes(no,Metric,color=names(three_dist_list[3])),size=size,alpha=alpha+0.3,linetype = "dashed")+
    
    
    
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)
}
#========================================================================================================================



#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc","mhd3","mhd3_collapsedCond")]
TF=( rowSums(df[,-c(6,7,8)])==5 )
sum(TF) #106



Original_pcc=cumsum( TF[order(df$pcc)])
Original_pcc_confusionMatrix_metrics=confusionMatrix_metrics(Original_pcc)


Ternary_mhd3=cumsum( TF[order(df$mhd3)]) 
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt with
Ternary_mhd3_confusionMatrix_metrics=confusionMatrix_metrics(Ternary_mhd3)


TernaryCondCollapsed_mhd3=cumsum( TF[order(df$mhd3_collapsedCond)])
##=>Some orders here are random because there are limited unique MHD distances. In the plot_metrics() this problem is properly dealt with
TernaryCondCollapsed_mhd3_confusionMatrix_metrics=confusionMatrix_metrics(TernaryCondCollapsed_mhd3)



three_dist_list=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
                     "Ternary data - MHD"=Ternary_mhd3_confusionMatrix_metrics,
                     "Ternary data with condition collapsed - MHD"=TernaryCondCollapsed_mhd3_confusionMatrix_metrics)



##Set no. of samples to be plotted
##samples=1:300
samples=1:5000
##samples=seq(from=0,to=7914231,by=1989); samples[1]=1

index_MHD3=df$mhd3 %>% table %>% cumsum 
index_collapsed_MHD3=df$mhd3_collapsedCond %>% table %>% cumsum

p1=plot_metrics(samples=samples,alpha=0.5,size=1,three_dist_list,
                index_MHD3,index_collapsed_MHD3,
                metric="sensitivity",y="Sensitivity") 
p2=plot_metrics(samples=samples,alpha=0.5,size=1,three_dist_list,
                index_MHD3,index_collapsed_MHD3,
                metric="specificity",y="specificity")
p3=plot_metrics(samples=samples,alpha=0.5,size=1,three_dist_list,
                index_MHD3,index_collapsed_MHD3,
                metric="precision",y="precision")
p4=plot_metrics(samples=samples,alpha=0.5,size=1,three_dist_list,
                index_MHD3,index_collapsed_MHD3,
                metric="accuracy",y="accuracy")


require(gridExtra)
p_all=grid.arrange(p1,p2,p3,p4, ncol=2) #When viewed in small window the dashed lines don't show properly

#dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
#ggsave(file=paste(dir_of_workingScript,"metrics_diffPhenotypeData.pdf",sep="/"),plot=p_all,device="pdf",width=3*6,height=2*6)

##When samples=1:300 or =1:5000 it cuts the edges so I changed the width and hieight a little bit
#ggsave(file=paste(dir_of_workingScript,"metrics_diffPhenotypeData_first5000.pdf",sep="/"),plot=p_all,device="pdf",width=3*7,height=2*7)
#ggsave(file=paste(dir_of_workingScript,"metrics_diffPhenotypeData_first300.pdf",sep="/"),plot=p_all,device="pdf",width=3*7,height=2*7)








