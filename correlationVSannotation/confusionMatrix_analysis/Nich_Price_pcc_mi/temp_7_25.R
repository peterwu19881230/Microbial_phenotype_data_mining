##This file is dependent on other files. Have to check

dist_list_Nich=list("Original data - PCC"=Original_pcc_confusionMatrix_metrics,
               "Ternary data - MI"=mi_ternary_confusionMatrix_metrics)

dist_list_Price=list("Price original data - PCC"=Price_pcc_confusionMatrix_metrics,
                     "Price Ternary data - MI"=Price_mi_ternary_confusionMatrix_metrics)




Price_samples=1:6211050

plot_PRC=function(samples,size,dist_list,legend_title="Type of data - Correlation",Cols){
  
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(recall=dist_list[[1]]$sensitivity[samples],precision=dist_list[[1]]$precision[samples]),aes(recall,precision,color=names(dist_list[1])),size=size,alpha=0.7)+
    geom_line(data=data.frame(recall=dist_list[[2]]$sensitivity[samples],precision=dist_list[[2]]$precision[samples]),aes(recall,precision,color=names(dist_list[2])),size=size,alpha=0.7)+
    
    
    #this is the negative control
    geom_abline(slope=0,intercept=(dist_list[[1]]$TP[1] + dist_list[[1]]$FN[1])/dim(dist_list[[1]])[1],linetype="dashed")+ #intercept= Positive/(Negative + Positive) = Positive/Total = (TP + FN)/total
    
    
    #scale_y_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    #scale_x_continuous(labels = function(x) paste0(round(x*100,3), "%"))+
    
    labs(x ="Recall (sensitivity)",y="Precision",aesthetic='custom text',color=legend_title)+
    #scale_x_continuous(expand = c(0,0),labels=comma)+
    scale_color_manual(values=Cols) 
}


prc_diffFitness=plot_PRC(Price_samples,size=1,dist_list_Nich,Cols=c("#09D38A","#5D46E0"))

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)

start.time=Sys.time()
ggsave(file=paste(dir_of_workingScript,"prc_diffFitness_Nich.png",sep="/"),plot=prc_diffFitness,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension
end.time=Sys.time()
end.time-start.time 


prc_diffFitness=plot_PRC(Price_samples,size=1,dist_list_Price,Cols=c("#832B05","#EC8D3A"))

start.time=Sys.time()
ggsave(file=paste(dir_of_workingScript,"prc_diffFitness_Price.png",sep="/"),plot=prc_diffFitness,width=9*1.5,height=4*1.5) #file type can be guessed by ggsave() by the extension
end.time=Sys.time()
end.time-start.time 





